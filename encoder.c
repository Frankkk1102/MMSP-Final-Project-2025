// encoder.c (C11) — MMSP Final Project JPEG-like pipeline (Mode 0..3)
// Build: gcc/clang -O2 -std=c11 -Wall -Wextra -pedantic -lm -o encoder encoder.c

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void die(const char *msg){ fprintf(stderr,"Error: %s\n", msg); exit(1); }
static int clampi(int v,int lo,int hi){ if(v<lo) return lo; if(v>hi) return hi; return v; }

// ---------- BMP 24-bit uncompressed ----------
#pragma pack(push,1)
typedef struct {
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} BFH;

typedef struct {
    uint32_t biSize;
    int32_t  biWidth;
    int32_t  biHeight;
    uint16_t biPlanes;
    uint16_t biBitCount;
    uint32_t biCompression;
    uint32_t biSizeImage;
    int32_t  biXPelsPerMeter;
    int32_t  biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
} BIH;
#pragma pack(pop)

typedef struct { int w,h; uint8_t *rgb; } ImageRGB;

static ImageRGB bmp_read_24(const char *path){
    FILE *fp=fopen(path,"rb"); if(!fp) die("cannot open bmp");
    BFH fh; BIH ih;
    if(fread(&fh,sizeof(fh),1,fp)!=1) die("bmp file header");
    if(fread(&ih,sizeof(ih),1,fp)!=1) die("bmp info header");
    if(fh.bfType!=0x4D42) die("not BMP");
    if(ih.biBitCount!=24 || ih.biCompression!=0) die("only 24-bit uncompressed BMP supported");

    int w=ih.biWidth;
    int h=(ih.biHeight>0)? ih.biHeight : -ih.biHeight;
    int bottom_up=(ih.biHeight>0);

    ImageRGB im; im.w=w; im.h=h;
    im.rgb=(uint8_t*)malloc((size_t)w*h*3); if(!im.rgb) die("OOM");

    int rowbytes=w*3;
    int pad=(4-(rowbytes%4))%4;

    fseek(fp, (long)fh.bfOffBits, SEEK_SET);
    for(int y=0;y<h;y++){
        int ry = bottom_up ? (h-1-y) : y;
        uint8_t *row = im.rgb + (size_t)ry*w*3;
        for(int x=0;x<w;x++){
            uint8_t bgr[3];
            if(fread(bgr,1,3,fp)!=3) die("bmp read pixel");
            row[3*x+0]=bgr[2]; // R
            row[3*x+1]=bgr[1]; // G
            row[3*x+2]=bgr[0]; // B
        }
        for(int p=0;p<pad;p++) (void)fgetc(fp);
    }
    fclose(fp);
    return im;
}

static void bmp_write_24(const char *path, const ImageRGB *im){
    FILE *fp=fopen(path,"wb"); if(!fp) die("cannot write bmp");

    int rowbytes=im->w*3;
    int pad=(4-(rowbytes%4))%4;
    uint32_t datasz=(uint32_t)((rowbytes+pad)*im->h);

    BFH fh={0};
    BIH ih={0};
    fh.bfType=0x4D42;
    fh.bfOffBits=(uint32_t)(sizeof(BFH)+sizeof(BIH));
    fh.bfSize=fh.bfOffBits+datasz;

    ih.biSize=40;
    ih.biWidth=im->w;
    ih.biHeight=im->h; // bottom-up
    ih.biPlanes=1;
    ih.biBitCount=24;
    ih.biCompression=0;
    ih.biSizeImage=datasz;

    fwrite(&fh,sizeof(fh),1,fp);
    fwrite(&ih,sizeof(ih),1,fp);

    for(int y=0;y<im->h;y++){
        int ry=im->h-1-y;
        const uint8_t *row = im->rgb + (size_t)ry*im->w*3;
        for(int x=0;x<im->w;x++){
            uint8_t bgr[3]={ row[3*x+2], row[3*x+1], row[3*x+0] };
            fwrite(bgr,1,3,fp);
        }
        for(int p=0;p<pad;p++) fputc(0,fp);
    }
    fclose(fp);
}

static void img_free(ImageRGB *im){ free(im->rgb); im->rgb=NULL; }

// ---------- Reversible Color Transform (lossless) ----------
typedef struct { int w,h; int16_t *Y,*Cb,*Cr; } ImageYCCi;

static ImageYCCi rgb_to_ycc_rct(const ImageRGB *im){
    ImageYCCi out; out.w=im->w; out.h=im->h;
    size_t n=(size_t)im->w*im->h;
    out.Y =(int16_t*)malloc(n*sizeof(int16_t));
    out.Cb=(int16_t*)malloc(n*sizeof(int16_t));
    out.Cr=(int16_t*)malloc(n*sizeof(int16_t));
    if(!out.Y||!out.Cb||!out.Cr) die("OOM");
    for(size_t i=0;i<n;i++){
        int R=im->rgb[3*i+0], G=im->rgb[3*i+1], B=im->rgb[3*i+2];
        int Y =(R + 2*G + B) >> 2;
        int Cb= B - G;
        int Cr= R - G;
        out.Y[i]=(int16_t)Y;
        out.Cb[i]=(int16_t)Cb;
        out.Cr[i]=(int16_t)Cr;
    }
    return out;
}

static void ycc_free(ImageYCCi *im){ free(im->Y); free(im->Cb); free(im->Cr); im->Y=im->Cb=im->Cr=NULL; }

// ---------- DCT/IDCT ----------
static void dct8x8(const double in[8][8], double out[8][8]){
    for(int u=0;u<8;u++){
        for(int v=0;v<8;v++){
            double sum=0.0;
            for(int x=0;x<8;x++){
                for(int y=0;y<8;y++){
                    sum += in[x][y] *
                           cos(((2*x+1)*u*M_PI)/16.0) *
                           cos(((2*y+1)*v*M_PI)/16.0);
                }
            }
            double cu=(u==0)?(1.0/sqrt(2.0)):1.0;
            double cv=(v==0)?(1.0/sqrt(2.0)):1.0;
            out[u][v]=0.25*cu*cv*sum;
        }
    }
}

static void idct8x8(const double in[8][8], double out[8][8]){
    for(int x=0;x<8;x++){
        for(int y=0;y<8;y++){
            double sum=0.0;
            for(int u=0;u<8;u++){
                for(int v=0;v<8;v++){
                    double cu=(u==0)?(1.0/sqrt(2.0)):1.0;
                    double cv=(v==0)?(1.0/sqrt(2.0)):1.0;
                    sum += cu*cv*in[u][v] *
                           cos(((2*x+1)*u*M_PI)/16.0) *
                           cos(((2*y+1)*v*M_PI)/16.0);
                }
            }
            out[x][y]=0.25*sum;
        }
    }
}

// ---------- Quantization tables (standard JPEG tables) ----------
static const int QT_Y[8][8]={
 {16,11,10,16,24,40,51,61},
 {12,12,14,19,26,58,60,55},
 {14,13,16,24,40,57,69,56},
 {14,17,22,29,51,87,80,62},
 {18,22,37,56,68,109,103,77},
 {24,35,55,64,81,104,113,92},
 {49,64,78,87,103,121,120,101},
 {72,92,95,98,112,100,103,99}
};

static const int QT_C[8][8]={
 {17,18,24,47,99,99,99,99},
 {18,21,26,66,99,99,99,99},
 {24,26,56,99,99,99,99,99},
 {47,66,99,99,99,99,99,99},
 {99,99,99,99,99,99,99,99},
 {99,99,99,99,99,99,99,99},
 {99,99,99,99,99,99,99,99},
 {99,99,99,99,99,99,99,99}
};

static void dump_qt_txt(const char *path, const int qt[8][8]){
    FILE *fp=fopen(path,"w"); if(!fp) die("open Qt txt");
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            fprintf(fp,"%d%c", qt[i][j], (j==7)?'\n':' ');
        }
    }
    fclose(fp);
}

// ---------- Block utilities ----------
static int16_t get_plane_pix(const int16_t *p,int w,int h,int x,int y){
    if(x<0) x=0; if(y<0) y=0;
    if(x>=w) x=w-1; if(y>=h) y=h-1;
    return p[(size_t)y*w + x];
}
static void get_block(const int16_t *p,int w,int h,int bx,int by,double out[8][8]){
    int x0=bx*8, y0=by*8;
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            out[i][j]=(double)get_plane_pix(p,w,h,x0+j,y0+i);
        }
    }
}

// ---------- ZigZag ----------
static const uint8_t ZZ[64]={
 0,1,8,16,9,2,3,10,
 17,24,32,25,18,11,4,5,
 12,19,26,33,40,48,41,34,
 27,20,13,6,7,14,21,28,
 35,42,49,56,57,50,43,36,
 29,22,15,23,30,37,44,51,
 58,59,52,45,38,31,39,46,
 53,60,61,54,47,55,62,63
};
static void block_to_zz(const int16_t b[8][8], int16_t out[64]){
    int16_t tmp[64];
    for(int i=0;i<8;i++) for(int j=0;j<8;j++) tmp[i*8+j]=b[i][j];
    for(int k=0;k<64;k++) out[k]=tmp[ZZ[k]];
}

// ---------- SQNR ----------
static double sqnr_db(double sig, double err){
    if(err<=0.0) return 1e9;
    return 10.0*log10(sig/err);
}

// ---------- RLE Pair ----------
typedef struct { uint8_t skip; int16_t val; } Pair;

// produce qblock (int16 quantized DCT) for a given plane/block
static void make_qblock(const int16_t *plane,int w,int h,int bx,int by,const int qt[8][8], int16_t q[8][8]){
    double blk[8][8], F[8][8];
    get_block(plane,w,h,bx,by,blk);
    dct8x8(blk,F);
    for(int u=0;u<8;u++){
        for(int v=0;v<8;v++){
            q[u][v]=(int16_t)llround(F[u][v]/(double)qt[u][v]);
        }
    }
}

// ---------- Mode 0 ----------
static void mode0(const char *bmp, const char *Rt, const char *Gt, const char *Bt, const char *dim){
    ImageRGB im=bmp_read_24(bmp);
    FILE *fr=fopen(Rt,"w"); if(!fr) die("R.txt");
    FILE *fg=fopen(Gt,"w"); if(!fg) die("G.txt");
    FILE *fb=fopen(Bt,"w"); if(!fb) die("B.txt");
    for(int y=0;y<im.h;y++){
        for(int x=0;x<im.w;x++){
            uint8_t R=im.rgb[(size_t)(y*im.w+x)*3+0];
            uint8_t G=im.rgb[(size_t)(y*im.w+x)*3+1];
            uint8_t B=im.rgb[(size_t)(y*im.w+x)*3+2];
            fprintf(fr,"%u%c",(unsigned)R,(x==im.w-1)?'\n':' ');
            fprintf(fg,"%u%c",(unsigned)G,(x==im.w-1)?'\n':' ');
            fprintf(fb,"%u%c",(unsigned)B,(x==im.w-1)?'\n':' ');
        }
    }
    fclose(fr); fclose(fg); fclose(fb);
    FILE *fd=fopen(dim,"w"); if(!fd) die("dim.txt");
    fprintf(fd,"%d %d\n", im.w, im.h);
    fclose(fd);
    img_free(&im);
}

// ---------- Mode 1 ----------
static void mode1(const char *bmp,
                  const char *QtYp,const char *QtCbp,const char *QtCrp,const char *dim,
                  const char *qFYp,const char *qFCbp,const char *qFCrp,
                  const char *eFYp,const char *eFCbp,const char *eFCrp)
{
    ImageRGB im=bmp_read_24(bmp);
    ImageYCCi ycc=rgb_to_ycc_rct(&im);

    dump_qt_txt(QtYp, QT_Y);
    dump_qt_txt(QtCbp, QT_C);
    dump_qt_txt(QtCrp, QT_C);

    FILE *fd=fopen(dim,"w"); if(!fd) die("dim.txt");
    fprintf(fd,"%d %d\n", im.w, im.h);
    fclose(fd);

    FILE *fqY=fopen(qFYp,"wb"); if(!fqY) die("qF_Y.raw");
    FILE *fqCb=fopen(qFCbp,"wb"); if(!fqCb) die("qF_Cb.raw");
    FILE *fqCr=fopen(qFCrp,"wb"); if(!fqCr) die("qF_Cr.raw");

    FILE *feY=fopen(eFYp,"wb"); if(!feY) die("eF_Y.raw");
    FILE *feCb=fopen(eFCbp,"wb"); if(!feCb) die("eF_Cb.raw");
    FILE *feCr=fopen(eFCrp,"wb"); if(!feCr) die("eF_Cr.raw");

    int nbx=(im.w+7)/8, nby=(im.h+7)/8;

    double sigY[64]={0}, errY[64]={0};
    double sigCb[64]={0}, errCb[64]={0};
    double sigCr[64]={0}, errCr[64]={0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const int16_t *plane = (comp==0)? ycc.Y : (comp==1)? ycc.Cb : ycc.Cr;
                const int (*qt)[8]   = (comp==0)? QT_Y : QT_C;

                double blk[8][8], F[8][8];
                get_block(plane, im.w, im.h, bx, by, blk);
                dct8x8(blk, F);

                int16_t q[8][8];
                float   e[8][8];

                for(int u=0;u<8;u++){
                    for(int v=0;v<8;v++){
                        double qv=(double)qt[u][v];
                        int16_t qij=(int16_t)llround(F[u][v]/qv);
                        double recon=(double)qij*qv;
                        double er=F[u][v]-recon;
                        q[u][v]=qij;
                        e[u][v]=(float)er;

                        int k=u*8+v;
                        if(comp==0){ sigY[k]+=F[u][v]*F[u][v]; errY[k]+=er*er; }
                        else if(comp==1){ sigCb[k]+=F[u][v]*F[u][v]; errCb[k]+=er*er; }
                        else { sigCr[k]+=F[u][v]*F[u][v]; errCr[k]+=er*er; }
                    }
                }

                FILE *fq=(comp==0)?fqY:(comp==1)?fqCb:fqCr;
                FILE *fe=(comp==0)?feY:(comp==1)?feCb:feCr;

                // write qF in row-major 64
                for(int u=0;u<8;u++) for(int v=0;v<8;v++) fwrite(&q[u][v], sizeof(int16_t), 1, fq);
                // write eF in same order
                for(int u=0;u<8;u++) for(int v=0;v<8;v++) fwrite(&e[u][v], sizeof(float), 1, fe);
            }
        }
    }

    fclose(fqY); fclose(fqCb); fclose(fqCr);
    fclose(feY); fclose(feCb); fclose(feCr);

    printf("SQNR_Y(64):\n");
    for(int k=0;k<64;k++) printf("%.6f%c", sqnr_db(sigY[k],errY[k]), (k%8==7)?'\n':' ');
    printf("SQNR_Cb(64):\n");
    for(int k=0;k<64;k++) printf("%.6f%c", sqnr_db(sigCb[k],errCb[k]), (k%8==7)?'\n':' ');
    printf("SQNR_Cr(64):\n");
    for(int k=0;k<64;k++) printf("%.6f%c", sqnr_db(sigCr[k],errCr[k]), (k%8==7)?'\n':' ');

    ycc_free(&ycc);
    img_free(&im);
}

// ---------- Mode 2: DPCM + ZigZag + RLE ----------
static void mode2_ascii(const char *bmp, const char *outtxt){
    ImageRGB im=bmp_read_24(bmp);
    ImageYCCi ycc=rgb_to_ycc_rct(&im);

    FILE *fp=fopen(outtxt,"w"); if(!fp) die("open rle_code.txt");
    fprintf(fp,"%d %d\n", im.w, im.h);

    int nbx=(im.w+7)/8, nby=(im.h+7)/8;
    int16_t prevDC[3]={0,0,0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const char *cname=(comp==0)?"Y":(comp==1)?"Cb":"Cr";
                const int16_t *plane = (comp==0)? ycc.Y : (comp==1)? ycc.Cb : ycc.Cr;
                const int (*qt)[8]   = (comp==0)? QT_Y : QT_C;

                int16_t q[8][8]; make_qblock(plane,im.w,im.h,bx,by,qt,q);
                int16_t zz[64]; block_to_zz(q,zz);

                // DPCM DC (store diff)
                int16_t dc=zz[0];
                int16_t diff=(int16_t)(dc - prevDC[comp]);
                prevDC[comp]=dc;
                zz[0]=diff;

                fprintf(fp,"(%d,%d, %s)", by, bx, cname);

                // RLE: (skip,value) for non-zeros, on the 64-length sequence
                int run=0;
                for(int k=0;k<64;k++){
                    if(zz[k]==0){
                        run++;
                    }else{
                        fprintf(fp," %d %d", run, (int)zz[k]);
                        run=0;
                    }
                }
                fprintf(fp,"\n");
            }
        }
    }

    fclose(fp);
    ycc_free(&ycc);
    img_free(&im);
}

static void mode2_binary(const char *bmp, const char *outbin){
    ImageRGB im=bmp_read_24(bmp);
    ImageYCCi ycc=rgb_to_ycc_rct(&im);

    FILE *fp=fopen(outbin,"wb"); if(!fp) die("open rle_code.bin");
    uint32_t w=(uint32_t)im.w, h=(uint32_t)im.h;
    uint32_t nbx=(uint32_t)((im.w+7)/8), nby=(uint32_t)((im.h+7)/8);
    fwrite(&w,4,1,fp); fwrite(&h,4,1,fp);
    fwrite(&nbx,4,1,fp); fwrite(&nby,4,1,fp);

    int16_t prevDC[3]={0,0,0};

    for(uint32_t by=0;by<nby;by++){
        for(uint32_t bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const int16_t *plane = (comp==0)? ycc.Y : (comp==1)? ycc.Cb : ycc.Cr;
                const int (*qt)[8]   = (comp==0)? QT_Y : QT_C;

                int16_t q[8][8]; make_qblock(plane,im.w,im.h,(int)bx,(int)by,qt,q);
                int16_t zz[64]; block_to_zz(q,zz);

                // DPCM DC
                int16_t dc=zz[0];
                int16_t diff=(int16_t)(dc - prevDC[comp]);
                prevDC[comp]=dc;
                zz[0]=diff;

                // Build pairs into temp
                Pair pairs[256];
                uint16_t np=0;
                int run=0;
                for(int k=0;k<64;k++){
                    if(zz[k]==0){
                        run++;
                    }else{
                        pairs[np].skip=(uint8_t)run;
                        pairs[np].val=zz[k];
                        np++;
                        run=0;
                    }
                }

                // write npairs then pairs
                fwrite(&np,2,1,fp);
                for(uint16_t i=0;i<np;i++){
                    fwrite(&pairs[i].skip,1,1,fp);
                    fwrite(&pairs[i].val,2,1,fp);
                }
            }
        }
    }

    fclose(fp);
    ycc_free(&ycc);
    img_free(&im);
}

// ---------- Huffman (Mode 3) ----------
typedef struct {
    uint32_t sym;     // packed (skip<<16)|(uint16)val
    uint32_t cnt;
    uint32_t code;    // bits in LSB
    uint8_t  len;
} Code;

typedef struct Node {
    uint32_t sym;
    uint32_t cnt;
    struct Node *l,*r;
} Node;

static Node* node_new(uint32_t sym,uint32_t cnt,Node* l,Node* r){
    Node *n=(Node*)calloc(1,sizeof(Node));
    if(!n) die("OOM");
    n->sym=sym; n->cnt=cnt; n->l=l; n->r=r;
    return n;
}
static void node_free(Node* n){
    if(!n) return;
    node_free(n->l);
    node_free(n->r);
    free(n);
}

static int cmp_u32(const void *a,const void *b){
    uint32_t x=*(const uint32_t*)a, y=*(const uint32_t*)b;
    return (x<y)?-1:(x>y);
}

static void build_codes_rec(Node *n, uint32_t code, uint8_t len, Code *out, int *idx){
    if(!n->l && !n->r){
        out[*idx].sym=n->sym;
        out[*idx].cnt=n->cnt;
        out[*idx].code=code;
        out[*idx].len=(len==0)?1:len; // single-symbol case
        (*idx)++;
        return;
    }
    if(n->l) build_codes_rec(n->l, (code<<1), (uint8_t)(len+1), out, idx);
    if(n->r) build_codes_rec(n->r, (code<<1)|1u, (uint8_t)(len+1), out, idx);
}

static int find_code_linear(const Code *codes,int nc,uint32_t sym,uint32_t *code,uint8_t *len){
    for(int i=0;i<nc;i++){
        if(codes[i].sym==sym){
            *code=codes[i].code;
            *len=codes[i].len;
            return 1;
        }
    }
    return 0;
}

static uint32_t pack_sym(uint8_t skip, int16_t val){
    return ((uint32_t)skip<<16) | (uint16_t)val;
}

static void bit_put(FILE *fp, uint32_t *buf, int *nbits, uint32_t code, uint8_t len){
    // append len bits of code (MSB-first) into buffer (left shift)
    *buf = (*buf<<len) | (code & ((len==32)?0xFFFFFFFFu:((1u<<len)-1u)));
    *nbits += (int)len;
    while(*nbits>=8){
        int shift=*nbits-8;
        uint8_t out=(uint8_t)((*buf>>shift)&0xFFu);
        fputc(out, fp);
        *nbits -= 8;
        *buf &= (shift==0)?0u:((1u<<shift)-1u);
    }
}
static void bit_flush(FILE *fp, uint32_t *buf, int *nbits){
    if(*nbits>0){
        uint8_t out=(uint8_t)((*buf<<(8-*nbits))&0xFFu);
        fputc(out, fp);
        *buf=0; *nbits=0;
    }
}

static long fsize(const char *path){
    FILE *fp=fopen(path,"rb"); if(!fp) return -1;
    fseek(fp,0,SEEK_END);
    long s=ftell(fp);
    fclose(fp);
    return s;
}

static void mode3(const char *bmp, const char *fmt, const char *codebook_txt, const char *outpath){
    // Build symbol stream from Mode2 pipeline but without writing rle file:
    ImageRGB im=bmp_read_24(bmp);
    ImageYCCi ycc=rgb_to_ycc_rct(&im);

    int nbx=(im.w+7)/8, nby=(im.h+7)/8;
    size_t blocks = (size_t)nbx*nby;

    // worst-case symbols: up to 64 nonzeros per block per channel => 64*blocks*3
    size_t cap = blocks*3*64;
    uint32_t *stream=(uint32_t*)malloc(cap*sizeof(uint32_t));
    if(!stream) die("OOM stream");
    size_t N=0;

    int16_t prevDC[3]={0,0,0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const int16_t *plane = (comp==0)? ycc.Y : (comp==1)? ycc.Cb : ycc.Cr;
                const int (*qt)[8]   = (comp==0)? QT_Y : QT_C;

                int16_t q[8][8]; make_qblock(plane,im.w,im.h,bx,by,qt,q);
                int16_t zz[64]; block_to_zz(q,zz);

                // DPCM DC
                int16_t dc=zz[0];
                int16_t diff=(int16_t)(dc - prevDC[comp]);
                prevDC[comp]=dc;
                zz[0]=diff;

                // RLE -> symbols
                int run=0;
                for(int k=0;k<64;k++){
                    if(zz[k]==0){
                        run++;
                    }else{
                        if(N>=cap) die("stream overflow (unexpected)");
                        stream[N++]=pack_sym((uint8_t)run, zz[k]);
                        run=0;
                    }
                }
                /* ---- ADD THIS: explicit End-Of-Block ---- */
                if(N >= cap) die("stream overflow (unexpected)");
                stream[N++] = pack_sym(255, 0);
                // note: no explicit EOB symbol; decoder stops by filling 64 slots rule
            }
        }
    }

    // Count frequencies by sorting a copy
    uint32_t *sorted=(uint32_t*)malloc(N*sizeof(uint32_t));
    if(!sorted) die("OOM sorted");
    memcpy(sorted, stream, N*sizeof(uint32_t));
    qsort(sorted, N, sizeof(uint32_t), cmp_u32);

    // unique symbols
    uint32_t *syms=(uint32_t*)malloc(N*sizeof(uint32_t));
    uint32_t *cnts=(uint32_t*)malloc(N*sizeof(uint32_t));
    if(!syms||!cnts) die("OOM freq");
    int ns=0;
    for(size_t i=0;i<N;){
        size_t j=i+1;
        while(j<N && sorted[j]==sorted[i]) j++;
        syms[ns]=sorted[i];
        cnts[ns]=(uint32_t)(j-i);
        ns++;
        i=j;
    }
    free(sorted);

    // Build Huffman tree (naive O(ns^2) merge, ok for homework)
    Node **nodes=(Node**)malloc((size_t)ns*sizeof(Node*));
    if(!nodes) die("OOM nodes");
    for(int i=0;i<ns;i++) nodes[i]=node_new(syms[i], cnts[i], NULL,NULL);

    int cur=ns;
    while(cur>1){
        int a=0,b=1;
        if(nodes[b]->cnt < nodes[a]->cnt){ int t=a; a=b; b=t; }
        for(int i=2;i<cur;i++){
            if(nodes[i]->cnt < nodes[a]->cnt){ b=a; a=i; }
            else if(nodes[i]->cnt < nodes[b]->cnt){ b=i; }
        }
        if(b<a){ int t=a; a=b; b=t; }
        Node *m=node_new(0, nodes[a]->cnt + nodes[b]->cnt, nodes[a], nodes[b]);
        nodes[a]=m;
        nodes[b]=nodes[cur-1];
        cur--;
    }
    Node *root=nodes[0];
    free(nodes);
    free(syms); free(cnts);

    // Generate code table
    Code *codes=(Code*)calloc((size_t)ns, sizeof(Code));
    if(!codes) die("OOM codes");
    int idx=0;
    build_codes_rec(root, 0, 0, codes, &idx);

    // Write codebook.txt (readable)
    FILE *fc=fopen(codebook_txt,"w"); if(!fc) die("open codebook.txt");
    fprintf(fc,"# columns: skip val count codeword\n");
    for(int i=0;i<idx;i++){
        uint8_t skip=(uint8_t)(codes[i].sym>>16);
        int16_t val=(int16_t)(uint16_t)(codes[i].sym & 0xFFFFu);
        fprintf(fc,"%u %d %u ", (unsigned)skip, (int)val, (unsigned)codes[i].cnt);
        for(int k=codes[i].len-1;k>=0;k--){
            fputc(((codes[i].code>>k)&1u)?'1':'0', fc);
        }
        fputc('\n', fc);
    }
    fclose(fc);

    if(strcmp(fmt,"ascii")==0){
        FILE *fo=fopen(outpath,"w"); if(!fo) die("open huffman_code.txt");
        fprintf(fo,"%d %d\n", im.w, im.h);
        // one symbol per line (bitstring)
        for(size_t i=0;i<N;i++){
            uint32_t code; uint8_t len;
            if(!find_code_linear(codes, idx, stream[i], &code, &len)) die("code lookup fail");
            for(int k=len-1;k>=0;k--) fputc(((code>>k)&1u)?'1':'0', fo);
            fputc('\n', fo);
        }
        fclose(fo);
    } else if(strcmp(fmt,"binary")==0){
        FILE *fo=fopen(outpath,"wb"); if(!fo) die("open huffman_code.bin");
        uint32_t w=(uint32_t)im.w, h=(uint32_t)im.h;
        fwrite(&w,4,1,fo); fwrite(&h,4,1,fo);

        uint64_t nbits=0;
        long pos=ftell(fo);
        fwrite(&nbits,8,1,fo); // placeholder

        uint32_t buf=0; int nbuf=0;
        for(size_t i=0;i<N;i++){
            uint32_t code; uint8_t len;
            if(!find_code_linear(codes, idx, stream[i], &code, &len)) die("code lookup fail");
            bit_put(fo,&buf,&nbuf,code,len);
            nbits += (uint64_t)len;
        }
        bit_flush(fo,&buf,&nbuf);

        long end=ftell(fo);
        fseek(fo,pos,SEEK_SET);
        fwrite(&nbits,8,1,fo);
        fseek(fo,end,SEEK_SET);
        fclose(fo);

        long orig=fsize(bmp), out=fsize(outpath);
        if(orig>0 && out>0){
            double ratio=(double)orig/(double)out;
            printf("Compression ratio (overall): %.6f (orig/out)\n", ratio);
            printf("Compression rate: %.2f%% (1 - out/orig)\n", 100.0*(1.0-(double)out/(double)orig));
        }
    } else {
        die("mode3: fmt must be ascii or binary");
    }

    free(stream);
    free(codes);
    node_free(root);
    ycc_free(&ycc);
    img_free(&im);
}

// ---------- main ----------
int main(int argc, char **argv){
    if(argc<2) die("usage: encoder <mode> ...");
    int mode=atoi(argv[1]);

    if(mode==0){
        if(argc!=7) die("encoder 0 in.bmp R.txt G.txt B.txt dim.txt");
        mode0(argv[2],argv[3],argv[4],argv[5],argv[6]);
        return 0;
    }
    if(mode==1){
        if(argc!=13) die("encoder 1 bmp Qt_Y Qt_Cb Qt_Cr dim qF_Y qF_Cb qF_Cr eF_Y eF_Cb eF_Cr");
        mode1(argv[2],argv[3],argv[4],argv[5],argv[6],
              argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
        return 0;
    }
    if(mode==2){
        if(argc!=5) die("encoder 2 bmp ascii|binary rle_code.(txt|bin)");
        if(strcmp(argv[3],"ascii")==0) mode2_ascii(argv[2],argv[4]);
        else if(strcmp(argv[3],"binary")==0) mode2_binary(argv[2],argv[4]);
        else die("mode2 fmt must be ascii or binary");
        return 0;
    }
    if(mode==3){
        if(argc!=6) die("encoder 3 bmp ascii|binary codebook.txt huffman_code.(txt|bin)");
        mode3(argv[2],argv[3],argv[4],argv[5]);
        return 0;
    }

    die("mode must be 0..3");
    return 0;
}
