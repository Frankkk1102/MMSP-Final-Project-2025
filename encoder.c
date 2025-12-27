// encoder.c — MMSP final JPEG pipeline (mode 0..3)
// Build: gcc -O2 -std=c11 -lm -o encoder encoder.c

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void die(const char *msg){ fprintf(stderr,"Error: %s\n", msg); exit(1); }

static long fsize(const char *path){
    FILE *fp=fopen(path,"rb"); if(!fp) return -1;
    fseek(fp,0,SEEK_END); long s=ftell(fp); fclose(fp); return s;
}

static int clampi(int v,int lo,int hi){ if(v<lo) return lo; if(v>hi) return hi; return v; }

// ---------------- BMP (24-bit, uncompressed) ----------------
#pragma pack(push,1)
typedef struct {
    uint16_t bfType;      // 'BM'
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
} BITMAPFILEHEADER;

typedef struct {
    uint32_t biSize;      // 40
    int32_t  biWidth;
    int32_t  biHeight;    // positive: bottom-up
    uint16_t biPlanes;    // 1
    uint16_t biBitCount;  // 24
    uint32_t biCompression; // 0
    uint32_t biSizeImage;
    int32_t  biXPelsPerMeter;
    int32_t  biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
} BITMAPINFOHEADER;
#pragma pack(pop)

typedef struct { int w,h; uint8_t *rgb; } ImageRGB;

static ImageRGB bmp_read_24(const char *path){
    FILE *fp=fopen(path,"rb"); if(!fp) die("cannot open bmp");
    BITMAPFILEHEADER fh; BITMAPINFOHEADER ih;
    if(fread(&fh,sizeof(fh),1,fp)!=1) die("bmp header");
    if(fread(&ih,sizeof(ih),1,fp)!=1) die("bmp info");
    if(fh.bfType!=0x4D42) die("not BMP");
    if(ih.biBitCount!=24 || ih.biCompression!=0) die("only 24-bit uncompressed BMP supported");
    int w=ih.biWidth;
    int h= (ih.biHeight>0)? ih.biHeight : -ih.biHeight;
    int bottom_up = (ih.biHeight>0);

    ImageRGB im; im.w=w; im.h=h;
    im.rgb=(uint8_t*)malloc((size_t)w*h*3); if(!im.rgb) die("OOM");

    int rowbytes = w*3;
    int pad = (4 - (rowbytes % 4)) % 4;

    fseek(fp, fh.bfOffBits, SEEK_SET);
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
        for(int p=0;p<pad;p++) fgetc(fp);
    }
    fclose(fp);
    return im;
}

static void img_free(ImageRGB *im){ free(im->rgb); im->rgb=NULL; }

// ---------------- Reversible Color Transform (RCT) ----------------
// Forward (lossless integer):
//   Y  = (R + 2G + B) >> 2
//   Cb = B - G
//   Cr = R - G
// Inverse:
//   G = Y - ((Cb + Cr) >> 2)
//   R = Cr + G
//   B = Cb + G
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

// ---------------- DCT / IDCT (double) ----------------
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

// ---------------- Quantization tables (fixed 8x8) ----------------
// You can replace these with your course's required tables if different.
static const int QT_Y_STD[8][8]={
 {16,11,10,16,24,40,51,61},
 {12,12,14,19,26,58,60,55},
 {14,13,16,24,40,57,69,56},
 {14,17,22,29,51,87,80,62},
 {18,22,37,56,68,109,103,77},
 {24,35,55,64,81,104,113,92},
 {49,64,78,87,103,121,120,101},
 {72,92,95,98,112,100,103,99}
};

static const int QT_C_STD[8][8]={
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
            fprintf(fp, "%d%c", qt[i][j], (j==7)?'\n':' ');
        }
    }
    fclose(fp);
}

// ---------------- ZigZag (JPEG order) ----------------
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

// ---------------- Block extraction (pad by edge) ----------------
static int16_t get_pix_plane(const int16_t *p,int w,int h,int x,int y){
    if(x<0) x=0; if(y<0) y=0;
    if(x>=w) x=w-1; if(y>=h) y=h-1;
    return p[(size_t)y*w + x];
}

static void get_block_i16(const int16_t *plane,int w,int h,int bx,int by,double out[8][8]){
    int x0=bx*8, y0=by*8;
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            int16_t v = get_pix_plane(plane,w,h,x0+j,y0+i);
            out[i][j]=(double)v;
        }
    }
}

// ---------------- SQNR helpers ----------------
static double sqnr_db(double sig_pow, double err_pow){
    if(err_pow<=0) return 1e9;
    return 10.0*log10(sig_pow/err_pow);
}

// ---------------- Mode 0 ----------------
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
            fprintf(fr,"%d%c",R,(x==im.w-1)?'\n':' ');
            fprintf(fg,"%d%c",G,(x==im.w-1)?'\n':' ');
            fprintf(fb,"%d%c",B,(x==im.w-1)?'\n':' ');
        }
    }

    fclose(fr); fclose(fg); fclose(fb);
    FILE *fd=fopen(dim,"w"); if(!fd) die("dim.txt");
    fprintf(fd,"%d %d\n", im.w, im.h);
    fclose(fd);

    img_free(&im);
}

// ---------------- Mode 1 (encoder) ----------------
static void mode1(const char *bmp,
                  const char *QtY, const char *QtCb, const char *QtCr, const char *dim,
                  const char *qFY, const char *qFCb, const char *qFCr,
                  const char *eFY, const char *eFCb, const char *eFCr)
{
    ImageRGB im=bmp_read_24(bmp);
    ImageYCCi ycc=rgb_to_ycc_rct(&im);

    dump_qt_txt(QtY, QT_Y_STD);
    dump_qt_txt(QtCb, QT_C_STD);
    dump_qt_txt(QtCr, QT_C_STD);

    FILE *fd=fopen(dim,"w"); if(!fd) die("dim.txt");
    fprintf(fd,"%d %d\n", im.w, im.h);
    fclose(fd);

    FILE *fqY=fopen(qFY,"wb"); if(!fqY) die("qF_Y.raw");
    FILE *fqCb=fopen(qFCb,"wb"); if(!fqCb) die("qF_Cb.raw");
    FILE *fqCr=fopen(qFCr,"wb"); if(!fqCr) die("qF_Cr.raw");

    FILE *feY=fopen(eFY,"wb"); if(!feY) die("eF_Y.raw");
    FILE *feCb=fopen(eFCb,"wb"); if(!feCb) die("eF_Cb.raw");
    FILE *feCr=fopen(eFCr,"wb"); if(!feCr) die("eF_Cr.raw");

    int nbx=(im.w+7)/8, nby=(im.h+7)/8;

    // accumulate SQNR per frequency (u,v) for each channel
    double sigY[64]={0}, errY[64]={0};
    double sigCb[64]={0}, errCb[64]={0};
    double sigCr[64]={0}, errCr[64]={0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const int16_t *plane = (comp==0)? ycc.Y : (comp==1)? ycc.Cb : ycc.Cr;
                const int (*qt)[8]   = (comp==0)? QT_Y_STD : QT_C_STD;

                double blk[8][8], F[8][8];
                get_block_i16(plane, im.w, im.h, bx, by, blk);

                // DCT
                dct8x8(blk, F);

                // Quant & error
                int16_t q[8][8];
                float   e[8][8];
                for(int u=0;u<8;u++){
                    for(int v=0;v<8;v++){
                        double qv = (double)qt[u][v];
                        int16_t qij = (int16_t)llround(F[u][v]/qv);
                        double recon = (double)qij * qv;
                        double err = F[u][v] - recon;
                        q[u][v]=qij;
                        e[u][v]=(float)err;

                        int k=u*8+v;
                        if(comp==0){ sigY[k]+=F[u][v]*F[u][v]; errY[k]+=err*err; }
                        else if(comp==1){ sigCb[k]+=F[u][v]*F[u][v]; errCb[k]+=err*err; }
                        else { sigCr[k]+=F[u][v]*F[u][v]; errCr[k]+=err*err; }
                    }
                }

                // write qF as int16 in row-major 64
                FILE *fq = (comp==0)?fqY:(comp==1)?fqCb:fqCr;
                FILE *fe = (comp==0)?feY:(comp==1)?feCb:feCr;

                for(int u=0;u<8;u++){
                    for(int v=0;v<8;v++){
                        fwrite(&q[u][v], sizeof(int16_t), 1, fq);
                    }
                }
                for(int u=0;u<8;u++){
                    for(int v=0;v<8;v++){
                        fwrite(&e[u][v], sizeof(float), 1, fe);
                    }
                }
            }
        }
    }

    fclose(fqY); fclose(fqCb); fclose(fqCr);
    fclose(feY); fclose(feCb); fclose(feCr);

    // Print 3x64 SQNR (dB)
    printf("SQNR_Y(64):\n");
    for(int k=0;k<64;k++){
        printf("%.6f%c", sqnr_db(sigY[k], errY[k]), (k%8==7)?'\n':' ');
    }
    printf("SQNR_Cb(64):\n");
    for(int k=0;k<64;k++){
        printf("%.6f%c", sqnr_db(sigCb[k], errCb[k]), (k%8==7)?'\n':' ');
    }
    printf("SQNR_Cr(64):\n");
    for(int k=0;k<64;k++){
        printf("%.6f%c", sqnr_db(sigCr[k], errCr[k]), (k%8==7)?'\n':' ');
    }

    ycc_free(&ycc);
    img_free(&im);
}

// ---------------- RLE stream (for mode2/3) ----------------
typedef struct { uint8_t skip; int16_t val; } Pair;

typedef struct {
    Pair *pairs;
    int npairs;
    int cap;
} PairList;

static void pl_init(PairList *pl){ pl->pairs=NULL; pl->npairs=0; pl->cap=0; }
static void pl_push(PairList *pl, uint8_t skip, int16_t val){
    if(pl->npairs==pl->cap){
        pl->cap = pl->cap? pl->cap*2 : 256;
        pl->pairs = (Pair*)realloc(pl->pairs, pl->cap*sizeof(Pair));
        if(!pl->pairs) die("OOM");
    }
    pl->pairs[pl->npairs].skip=skip;
    pl->pairs[pl->npairs].val=val;
    pl->npairs++;
}
static void pl_free(PairList *pl){ free(pl->pairs); pl->pairs=NULL; pl->npairs=pl->cap=0; }

// build qF (int16[8][8]) for each block&channel then apply DPCM+ZigZag+RLE on 64 sequence
static void build_rle_all(const char *bmp, PairList *outY, PairList *outCb, PairList *outCr, int *W, int *H,
                          int verbose_ascii, FILE *fascii)
{
    ImageRGB im=bmp_read_24(bmp);
    ImageYCCi ycc=rgb_to_ycc_rct(&im);
    *W=im.w; *H=im.h;

    int nbx=(im.w+7)/8, nby=(im.h+7)/8;

    int16_t prevDC_Y=0, prevDC_Cb=0, prevDC_Cr=0;

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const int16_t *plane = (comp==0)? ycc.Y : (comp==1)? ycc.Cb : ycc.Cr;
                const int (*qt)[8]   = (comp==0)? QT_Y_STD : QT_C_STD;

                double blk[8][8], F[8][8];
                get_block_i16(plane, im.w, im.h, bx, by, blk);
                dct8x8(blk, F);

                int16_t q[8][8];
                for(int u=0;u<8;u++){
                    for(int v=0;v<8;v++){
                        q[u][v]=(int16_t)llround(F[u][v]/(double)qt[u][v]);
                    }
                }

                int16_t zz[64];
                block_to_zz(q, zz);

                // DPCM on DC: replace zz[0] with diff
                int16_t *pPrev = (comp==0)?&prevDC_Y:(comp==1)?&prevDC_Cb:&prevDC_Cr;
                int16_t dc = zz[0];
                zz[0] = (int16_t)(dc - *pPrev);
                *pPrev = dc;

                // RLE on full 64 (including DC diff at index0)
                PairList *pl = (comp==0)?outY:(comp==1)?outCb:outCr;

                // For ASCII output, print header per block/channel then pairs
                if(verbose_ascii){
                    const char *cname = (comp==0)?"Y":(comp==1)?"Cb":"Cr";
                    fprintf(fascii,"(%d,%d, %s)", by, bx, cname);
                }

                int run=0;
                for(int k=0;k<64;k++){
                    if(zz[k]==0){
                        run++;
                    }else{
                        uint8_t skip=(uint8_t)run;
                        pl_push(pl, skip, zz[k]);
                        if(verbose_ascii) fprintf(fascii," %u %d", skip, (int)zz[k]);
                        run=0;
                    }
                }
                if(verbose_ascii) fprintf(fascii,"\n");
            }
        }
    }

    ycc_free(&ycc);
    img_free(&im);
}

// ---------------- Mode 2 (encoder) ----------------
static void mode2(const char *bmp, const char *fmt, const char *outpath){
    PairList y,cb,cr; pl_init(&y); pl_init(&cb); pl_init(&cr);
    int W,H;

    if(strcmp(fmt,"ascii")==0){
        FILE *fa=fopen(outpath,"w"); if(!fa) die("open rle_code.txt");
        build_rle_all(bmp,&y,&cb,&cr,&W,&H,1,fa);
        // first row size should be first, per spec. We wrote later. Fix: prepend by rewriting.
        // easiest: create temp then rewrite. For simplicity: do it right now: use temp.
        fclose(fa);

        // Re-generate with correct first row:
        // We'll write size first, then regenerate (without storing extra again).
        // (cost is fine for homework)
        pl_free(&y); pl_free(&cb); pl_free(&cr);
        pl_init(&y); pl_init(&cb); pl_init(&cr);

        fa=fopen(outpath,"w"); if(!fa) die("open rle_code.txt");
        // Need W,H first; read from bmp quickly:
        ImageRGB im=bmp_read_24(bmp); W=im.w; H=im.h; img_free(&im);
        fprintf(fa,"%d %d\n", W, H);
        build_rle_all(bmp,&y,&cb,&cr,&W,&H,1,fa);
        fclose(fa);
    } else if(strcmp(fmt,"binary")==0){
        build_rle_all(bmp,&y,&cb,&cr,&W,&H,0,NULL);

        FILE *fb=fopen(outpath,"wb"); if(!fb) die("open rle_code.bin");
        uint32_t w=(uint32_t)W,h=(uint32_t)H;
        fwrite(&w,4,1,fb); fwrite(&h,4,1,fb);

        // For binary, we must write per-block/channel, but we stored as one flat list (still decodable if decoder mirrors).
        // We'll store three lists in order: Y then Cb then Cr with their lengths.
        // (You said format self-defined, decoder can understand.)
        uint32_t ny=(uint32_t)y.npairs, ncb=(uint32_t)cb.npairs, ncr=(uint32_t)cr.npairs;
        fwrite(&ny,4,1,fb); fwrite(&ncb,4,1,fb); fwrite(&ncr,4,1,fb);

        fwrite(y.pairs, sizeof(Pair), y.npairs, fb);
        fwrite(cb.pairs,sizeof(Pair), cb.npairs, fb);
        fwrite(cr.pairs,sizeof(Pair), cr.npairs, fb);
        fclose(fb);

        // Print compression ratio
        long orig=fsize(bmp);
        long out =fsize(outpath);
        if(orig>0 && out>0){
            double ratio = (double)orig/(double)out;
            printf("Compression ratio (overall): %.6f (orig/out)\n", ratio);
            printf("Compression rate: %.2f%% (1 - out/orig)\n", 100.0*(1.0-(double)out/(double)orig));
        }
    } else {
        die("mode2: fmt must be ascii or binary");
    }

    pl_free(&y); pl_free(&cb); pl_free(&cr);
}

// ---------------- Huffman (Mode 3) ----------------
typedef struct { uint32_t sym; uint32_t cnt; } SymCount;

typedef struct Node {
    uint32_t sym;
    uint32_t cnt;
    struct Node *l,*r;
} Node;

static Node* node_new(uint32_t sym,uint32_t cnt,Node* l,Node* r){
    Node* n=(Node*)calloc(1,sizeof(Node));
    if(!n) die("OOM");
    n->sym=sym; n->cnt=cnt; n->l=l; n->r=r;
    return n;
}
static void node_free(Node* n){ if(!n) return; node_free(n->l); node_free(n->r); free(n); }

static int cmp_sc(const void*a,const void*b){
    const SymCount* x=(const SymCount*)a;
    const SymCount* y=(const SymCount*)b;
    if(x->cnt<y->cnt) return -1;
    if(x->cnt>y->cnt) return  1;
    return (x->sym<y->sym)?-1:(x->sym>y->sym);
}

typedef struct { uint32_t sym; uint32_t code; uint8_t len; uint32_t cnt; } Code;

static void build_codes_rec(Node* n, uint32_t code, uint8_t len, Code *out, int *idx){
    if(!n->l && !n->r){
        out[*idx].sym=n->sym;
        out[*idx].code=code;
        out[*idx].len=len?len:1; // avoid zero-length if only one symbol
        out[*idx].cnt=n->cnt;
        (*idx)++;
        return;
    }
    if(n->l) build_codes_rec(n->l, (code<<1), (uint8_t)(len+1), out, idx);
    if(n->r) build_codes_rec(n->r, (code<<1)|1u, (uint8_t)(len+1), out, idx);
}

static void bit_put(FILE *fp, uint32_t *buf, int *nbits, uint32_t code, uint8_t len){
    *buf = (*buf<<len) | (code & ((len==32)?0xFFFFFFFFu:((1u<<len)-1u)));
    *nbits += len;
    while(*nbits>=8){
        int shift=*nbits-8;
        uint8_t out=(uint8_t)((*buf>>shift)&0xFFu);
        fputc(out, fp);
        *nbits-=8;
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

// symbol definition for Huffman:
// pack (skip,val) into 32-bit: sym = (skip<<16) | (uint16_t)val
static uint32_t pack_sym(uint8_t skip, int16_t val){
    return ((uint32_t)skip<<16) | (uint16_t)val;
}

static void mode3(const char *bmp, const char *fmt, const char *codebook_txt, const char *outpath){
    PairList y,cb,cr; pl_init(&y); pl_init(&cb); pl_init(&cr);
    int W,H;
    build_rle_all(bmp,&y,&cb,&cr,&W,&H,0,NULL);

    // Build symbol stream: concat Y then Cb then Cr (decoder will do same)
    size_t N = (size_t)y.npairs + cb.npairs + cr.npairs;
    uint32_t *stream = (uint32_t*)malloc(N*sizeof(uint32_t));
    if(!stream) die("OOM");
    size_t t=0;
    for(int i=0;i<y.npairs;i++)  stream[t++]=pack_sym(y.pairs[i].skip,  y.pairs[i].val);
    for(int i=0;i<cb.npairs;i++) stream[t++]=pack_sym(cb.pairs[i].skip, cb.pairs[i].val);
    for(int i=0;i<cr.npairs;i++) stream[t++]=pack_sym(cr.pairs[i].skip, cr.pairs[i].val);

    // Count frequencies
    // Simple approach: sort stream then run-length count
    uint32_t *sorted=(uint32_t*)malloc(N*sizeof(uint32_t));
    if(!sorted) die("OOM");
    memcpy(sorted, stream, N*sizeof(uint32_t));
    // qsort needs comparator
    int cmp_u32(const void*a,const void*b){
        uint32_t x=*(const uint32_t*)a, y=*(const uint32_t*)b;
        return (x<y)?-1:(x>y);
    }
    qsort(sorted,N,sizeof(uint32_t),cmp_u32);

    SymCount *sc=(SymCount*)malloc(N*sizeof(SymCount));
    if(!sc) die("OOM");
    int ns=0;
    for(size_t i=0;i<N;){
        size_t j=i+1;
        while(j<N && sorted[j]==sorted[i]) j++;
        sc[ns].sym=sorted[i];
        sc[ns].cnt=(uint32_t)(j-i);
        ns++;
        i=j;
    }
    free(sorted);

    // Build Huffman tree using naive O(ns^2) merging (fine for homework)
    // Use an array of Node* as a multiset sorted by cnt
    Node **nodes=(Node**)malloc(ns*sizeof(Node*));
    if(!nodes) die("OOM");
    for(int i=0;i<ns;i++) nodes[i]=node_new(sc[i].sym, sc[i].cnt, NULL,NULL);

    int cur=ns;
    while(cur>1){
        // find two smallest
        int a=0,b=1;
        if(nodes[b]->cnt < nodes[a]->cnt){ int tmp=a;a=b;b=tmp; }
        for(int i=2;i<cur;i++){
            if(nodes[i]->cnt < nodes[a]->cnt){ b=a; a=i; }
            else if(nodes[i]->cnt < nodes[b]->cnt){ b=i; }
        }
        if(b<a){ int tmp=a;a=b;b=tmp; }

        Node *m = node_new(0, nodes[a]->cnt + nodes[b]->cnt, nodes[a], nodes[b]);

        // remove b then a, insert m
        nodes[a]=m;
        nodes[b]=nodes[cur-1];
        cur--;
    }
    Node *root = nodes[0];
    free(nodes);

    // Generate codes
    Code *codes=(Code*)calloc(ns,sizeof(Code));
    if(!codes) die("OOM");
    int idx=0;
    build_codes_rec(root, 0, 0, codes, &idx);

    // Write codebook.txt (readable)
    FILE *fc=fopen(codebook_txt,"w"); if(!fc) die("codebook.txt");
    fprintf(fc,"# symbol format: sym=(skip<<16)|uint16(val)\n");
    fprintf(fc,"# columns: skip val count codeword\n");
    for(int i=0;i<idx;i++){
        uint8_t skip = (uint8_t)(codes[i].sym>>16);
        int16_t val  = (int16_t)(uint16_t)(codes[i].sym & 0xFFFFu);
        fprintf(fc,"%u %d %u ", skip, (int)val, (unsigned)codes[i].cnt);
        for(int k=codes[i].len-1;k>=0;k--){
            fputc(((codes[i].code>>k)&1u)?'1':'0', fc);
        }
        fputc('\n', fc);
    }
    fclose(fc);

    // Helper to lookup code by symbol (linear; okay for homework)
    auto int find_code(uint32_t sym, uint32_t *code, uint8_t *len){
        for(int i=0;i<idx;i++){
            if(codes[i].sym==sym){ *code=codes[i].code; *len=codes[i].len; return 1; }
        }
        return 0;
    }

    if(strcmp(fmt,"ascii")==0){
        FILE *fo=fopen(outpath,"w"); if(!fo) die("huffman_code.txt");
        fprintf(fo,"%d %d\n", W, H);
        // dump bitstream as '0''1' chars (readable)
        for(size_t i=0;i<N;i++){
            uint32_t code; uint8_t len;
            if(!find_code(stream[i], &code, &len)) die("code lookup fail");
            for(int k=len-1;k>=0;k--) fputc(((code>>k)&1u)?'1':'0', fo);
            fputc('\n', fo); // one symbol per line (easiest for decoder)
        }
        fclose(fo);
    } else if(strcmp(fmt,"binary")==0){
        FILE *fo=fopen(outpath,"wb"); if(!fo) die("huffman_code.bin");
        uint32_t w=(uint32_t)W,h=(uint32_t)H;
        fwrite(&w,4,1,fo); fwrite(&h,4,1,fo);

        // write nbits placeholder then bitstream
        uint64_t nbits=0;
        long pos_bits = ftell(fo);
        fwrite(&nbits,8,1,fo);

        uint32_t buf=0; int nbuf=0;
        for(size_t i=0;i<N;i++){
            uint32_t code; uint8_t len;
            if(!find_code(stream[i], &code, &len)) die("code lookup fail");
            bit_put(fo,&buf,&nbuf,code,len);
            nbits += len;
        }
        bit_flush(fo,&buf,&nbuf);

        // patch nbits
        long end = ftell(fo);
        fseek(fo,pos_bits,SEEK_SET);
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
    free(sc);
    free(codes);
    node_free(root);
    pl_free(&y); pl_free(&cb); pl_free(&cr);
}

int main(int argc, char **argv){
    if(argc<2) die("usage: encoder <mode> ...");

    int mode=atoi(argv[1]);

    if(mode==0){
        if(argc!=7) die("encoder 0 in.bmp R.txt G.txt B.txt dim.txt");
        mode0(argv[2],argv[3],argv[4],argv[5],argv[6]);
        return 0;
    }
    if(mode==1){
        if(argc!=13) die("encoder 1 bmp QtY QtCb QtCr dim qFY qFCb qFCr eFY eFCb eFCr");
        mode1(argv[2],argv[3],argv[4],argv[5],argv[6],
              argv[7],argv[8],argv[9],argv[10],argv[11],argv[12]);
        return 0;
    }
    if(mode==2){
        if(argc!=5) die("encoder 2 bmp ascii|binary out");
        mode2(argv[2],argv[3],argv[4]);
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
