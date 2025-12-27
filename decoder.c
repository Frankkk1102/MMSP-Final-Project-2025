// decoder.c (C11) — MMSP Final Project JPEG-like pipeline (Mode 0..3)
// Build: gcc/clang -O2 -std=c11 -Wall -Wextra -pedantic -lm -o decoder decoder.c

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
            row[3*x+0]=bgr[2];
            row[3*x+1]=bgr[1];
            row[3*x+2]=bgr[0];
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
    ih.biSizeImage=0;

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

// ---------- Reversible Color Transform inverse ----------
typedef struct { int w,h; int16_t *Y,*Cb,*Cr; } ImageYCCi;

static ImageRGB ycc_to_rgb_rct(const ImageYCCi *ycc){
    ImageRGB out; out.w=ycc->w; out.h=ycc->h;
    size_t n=(size_t)out.w*out.h;
    out.rgb=(uint8_t*)malloc(n*3); if(!out.rgb) die("OOM");
    for(size_t i=0;i<n;i++){
        int Y = ycc->Y[i];
        int Cb= ycc->Cb[i];
        int Cr= ycc->Cr[i];
        int G = Y - ((Cb + Cr) >> 2);
        int R = Cr + G;
        int B = Cb + G;
        out.rgb[3*i+0]=(uint8_t)clampi(R,0,255);
        out.rgb[3*i+1]=(uint8_t)clampi(G,0,255);
        out.rgb[3*i+2]=(uint8_t)clampi(B,0,255);
    }
    return out;
}

static void ycc_free(ImageYCCi *im){ free(im->Y); free(im->Cb); free(im->Cr); im->Y=im->Cb=im->Cr=NULL; }

// ---------- DCT/IDCT ----------
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

// ---------- QT read ----------
static void read_qt_txt(const char *path, int qt[8][8]){
    FILE *fp=fopen(path,"r"); if(!fp) die("open Qt txt");
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            if(fscanf(fp,"%d",&qt[i][j])!=1) die("Qt parse");
        }
    }
    fclose(fp);
}

// ---------- ZigZag inverse ----------
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
static void zz_to_block(const int16_t in[64], int16_t b[8][8]){
    int16_t tmp[64]={0};
    for(int k=0;k<64;k++) tmp[ZZ[k]]=in[k];
    for(int i=0;i<8;i++) for(int j=0;j<8;j++) b[i][j]=tmp[i*8+j];
}

// ---------- Put block into plane ----------
static void put_block_i16(int16_t *plane,int w,int h,int bx,int by,const double blk[8][8]){
    int x0=bx*8, y0=by*8;
    for(int i=0;i<8;i++){
        for(int j=0;j<8;j++){
            int x=x0+j, y=y0+i;
            if(x<w && y<h){
                plane[(size_t)y*w + x] = (int16_t)llround(blk[i][j]);
            }
        }
    }
}

// ---------- Mode 0 ----------
static void mode0(const char *outbmp, const char *Rt, const char *Gt, const char *Bt, const char *dim){
    int W,H;
    FILE *fd=fopen(dim,"r"); if(!fd) die("dim.txt");
    if(fscanf(fd,"%d %d",&W,&H)!=2) die("dim parse");
    fclose(fd);

    ImageRGB im; im.w=W; im.h=H;
    im.rgb=(uint8_t*)malloc((size_t)W*H*3); if(!im.rgb) die("OOM");

    FILE *fr=fopen(Rt,"r"); if(!fr) die("R.txt");
    FILE *fg=fopen(Gt,"r"); if(!fg) die("G.txt");
    FILE *fb=fopen(Bt,"r"); if(!fb) die("B.txt");

    for(int y=0;y<H;y++){
        for(int x=0;x<W;x++){
            int r,g,b;
            if(fscanf(fr,"%d",&r)!=1) die("R parse");
            if(fscanf(fg,"%d",&g)!=1) die("G parse");
            if(fscanf(fb,"%d",&b)!=1) die("B parse");
            im.rgb[(size_t)(y*W+x)*3+0]=(uint8_t)clampi(r,0,255);
            im.rgb[(size_t)(y*W+x)*3+1]=(uint8_t)clampi(g,0,255);
            im.rgb[(size_t)(y*W+x)*3+2]=(uint8_t)clampi(b,0,255);
        }
    }

    fclose(fr); fclose(fg); fclose(fb);
    bmp_write_24(outbmp,&im);
    img_free(&im);
}

// ---------- Pixel SQNR (RGB) ----------
static void print_pixel_sqnr_rgb(const ImageRGB *orig, const ImageRGB *rec){
    double sig[3]={0}, err[3]={0};
    size_t n=(size_t)orig->w*orig->h;
    for(size_t i=0;i<n;i++){
        for(int c=0;c<3;c++){
            double s=(double)orig->rgb[3*i+c];
            double d=(double)orig->rgb[3*i+c] - (double)rec->rgb[3*i+c];
            sig[c]+=s*s;
            err[c]+=d*d;
        }
    }
    for(int c=0;c<3;c++){
        double sq = (err[c]<=0)?1e9:10.0*log10(sig[c]/err[c]);
        const char *name=(c==0)?"R":(c==1)?"G":"B";
        printf("Pixel SQNR %s: %.6f dB\n", name, sq);
    }
}

// ---------- Mode 1 decode from qF (+ optional eF) ----------
static void decode_mode1(const char *outbmp,
                         int W,int H,
                         const int qtY[8][8], const int qtC[8][8],
                         FILE *fqY, FILE *fqCb, FILE *fqCr,
                         FILE *feY, FILE *feCb, FILE *feCr,
                         int use_error)
{
    ImageYCCi ycc; ycc.w=W; ycc.h=H;
    size_t n=(size_t)W*H;
    ycc.Y =(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cb=(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cr=(int16_t*)calloc(n,sizeof(int16_t));
    if(!ycc.Y||!ycc.Cb||!ycc.Cr) die("OOM");

    int nbx=(W+7)/8, nby=(H+7)/8;

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                const int (*qt)[8]=(comp==0)?qtY:qtC;
                FILE *fq=(comp==0)?fqY:(comp==1)?fqCb:fqCr;
                FILE *fe=(comp==0)?feY:(comp==1)?feCb:feCr;

                int16_t q[8][8];
                float   e[8][8];

                for(int u=0;u<8;u++) for(int v=0;v<8;v++){
                    if(fread(&q[u][v],sizeof(int16_t),1,fq)!=1) die("qF read");
                }
                if(use_error){
                    for(int u=0;u<8;u++) for(int v=0;v<8;v++){
                        if(fread(&e[u][v],sizeof(float),1,fe)!=1) die("eF read");
                    }
                }

                double F[8][8];
                for(int u=0;u<8;u++){
                    for(int v=0;v<8;v++){
                        double base=(double)q[u][v]*(double)qt[u][v];
                        F[u][v]= use_error ? (base + (double)e[u][v]) : base;
                    }
                }

                double blk[8][8];
                idct8x8(F,blk);

                int16_t *plane=(comp==0)?ycc.Y:(comp==1)?ycc.Cb:ycc.Cr;
                put_block_i16(plane,W,H,bx,by,blk);
            }
        }
    }

    ImageRGB rgb=ycc_to_rgb_rct(&ycc);
    bmp_write_24(outbmp,&rgb);
    img_free(&rgb);
    ycc_free(&ycc);
}

// ---------- Mode 2 & 3 shared: RLE pair ----------
typedef struct { uint8_t skip; int16_t val; } Pair;

static void rle_pairs_to_zz(int16_t zz[64], const Pair *pairs, int npairs){
    memset(zz,0,64*sizeof(int16_t));
    int k=0;
    for(int i=0;i<npairs;i++){
        k += (int)pairs[i].skip;
        if(k>=64) break;
        zz[k]=pairs[i].val;
        k++;
        if(k>=64) break;
    }
}

static void inv_dpcm_dc(int16_t zz[64], int16_t *prevDC){
    int16_t diff=zz[0];
    int16_t dc=(int16_t)(diff + *prevDC);
    zz[0]=dc;
    *prevDC=dc;
}

static void reconstruct_from_rle_blocks(const char *outbmp, int W,int H,
                                       const int qtY[8][8], const int qtC[8][8],
                                       int nbx,int nby,
                                       FILE *get_next_block_pairs, int is_dummy) {
    (void)outbmp; (void)W; (void)H; (void)qtY; (void)qtC; (void)nbx; (void)nby; (void)get_next_block_pairs; (void)is_dummy;
    // placeholder (not used)
}

// ---------- Mode 2 decode ASCII ----------
static void mode2_ascii(const char *outbmp, const char *txtpath){
    FILE *fp=fopen(txtpath,"r"); if(!fp) die("open rle_code.txt");
    int W,H;
    if(fscanf(fp,"%d %d",&W,&H)!=2) die("size row");
    int nbx=(W+7)/8, nby=(H+7)/8;

    // QT must match encoder’s fixed QT (mode2 spec doesn’t pass Qt files)
    static const int qtY[8][8]={
     {16,11,10,16,24,40,51,61},{12,12,14,19,26,58,60,55},{14,13,16,24,40,57,69,56},{14,17,22,29,51,87,80,62},
     {18,22,37,56,68,109,103,77},{24,35,55,64,81,104,113,92},{49,64,78,87,103,121,120,101},{72,92,95,98,112,100,103,99}
    };
    static const int qtC[8][8]={
     {17,18,24,47,99,99,99,99},{18,21,26,66,99,99,99,99},{24,26,56,99,99,99,99,99},{47,66,99,99,99,99,99,99},
     {99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99}
    };

    // consume rest of first line
    char line[1<<16];
    fgets(line,sizeof(line),fp);

    ImageYCCi ycc; ycc.w=W; ycc.h=H;
    size_t n=(size_t)W*H;
    ycc.Y =(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cb=(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cr=(int16_t*)calloc(n,sizeof(int16_t));
    if(!ycc.Y||!ycc.Cb||!ycc.Cr) die("OOM");

    int16_t prevDC[3]={0,0,0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                if(!fgets(line,sizeof(line),fp)) die("EOF rle_code.txt");
                // format: (m,n, X) skip val skip val ...
                char *p=strchr(line,')');
                if(!p) die("bad rle line");
                p++;

                Pair pairs[256]; int np=0;
                while(1){
                    while(*p==' '||*p=='\t') p++;
                    if(*p=='\0'||*p=='\n'||*p=='\r') break;
                    int skip,val;
                    if(sscanf(p,"%d %d",&skip,&val)!=2) break;
                    pairs[np].skip=(uint8_t)skip;
                    pairs[np].val=(int16_t)val;
                    np++;

                    // advance p twice
                    while(*p && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r') p++;
                    while(*p==' '||*p=='\t') p++;
                    while(*p && *p!=' ' && *p!='\t' && *p!='\n' && *p!='\r') p++;
                }

                int16_t zz[64]; rle_pairs_to_zz(zz,pairs,np);
                inv_dpcm_dc(zz, &prevDC[comp]);

                int16_t qblk[8][8]; zz_to_block(zz,qblk);

                const int (*qt)[8]=(comp==0)?qtY:qtC;
                double F[8][8];
                for(int u=0;u<8;u++) for(int v=0;v<8;v++) F[u][v]=(double)qblk[u][v]*(double)qt[u][v];

                double blk[8][8];
                idct8x8(F,blk);

                int16_t *plane=(comp==0)?ycc.Y:(comp==1)?ycc.Cb:ycc.Cr;
                put_block_i16(plane,W,H,bx,by,blk);
            }
        }
    }

    fclose(fp);
    ImageRGB rgb=ycc_to_rgb_rct(&ycc);
    bmp_write_24(outbmp,&rgb);
    img_free(&rgb);
    ycc_free(&ycc);
}

// ---------- Mode 2 decode Binary (our format) ----------
static void mode2_binary(const char *outbmp, const char *binpath){
    FILE *fp=fopen(binpath,"rb"); if(!fp) die("open rle_code.bin");
    uint32_t W,H,nbx,nby;
    if(fread(&W,4,1,fp)!=1) die("bin read W");
    if(fread(&H,4,1,fp)!=1) die("bin read H");
    if(fread(&nbx,4,1,fp)!=1) die("bin read nbx");
    if(fread(&nby,4,1,fp)!=1) die("bin read nby");

    static const int qtY[8][8]={
     {16,11,10,16,24,40,51,61},{12,12,14,19,26,58,60,55},{14,13,16,24,40,57,69,56},{14,17,22,29,51,87,80,62},
     {18,22,37,56,68,109,103,77},{24,35,55,64,81,104,113,92},{49,64,78,87,103,121,120,101},{72,92,95,98,112,100,103,99}
    };
    static const int qtC[8][8]={
     {17,18,24,47,99,99,99,99},{18,21,26,66,99,99,99,99},{24,26,56,99,99,99,99,99},{47,66,99,99,99,99,99,99},
     {99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99}
    };

    ImageYCCi ycc; ycc.w=(int)W; ycc.h=(int)H;
    size_t n=(size_t)W*H;
    ycc.Y =(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cb=(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cr=(int16_t*)calloc(n,sizeof(int16_t));
    if(!ycc.Y||!ycc.Cb||!ycc.Cr) die("OOM");

    int16_t prevDC[3]={0,0,0};

    for(uint32_t by=0;by<nby;by++){
        for(uint32_t bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                uint16_t np;
                if(fread(&np,2,1,fp)!=1) die("bin read npairs");
                Pair pairs[256]; if(np>256) die("npairs too big");
                for(uint16_t i=0;i<np;i++){
                    if(fread(&pairs[i].skip,1,1,fp)!=1) die("bin read skip");
                    if(fread(&pairs[i].val,2,1,fp)!=1) die("bin read val");
                }

                int16_t zz[64]; rle_pairs_to_zz(zz,pairs,(int)np);
                inv_dpcm_dc(zz, &prevDC[comp]);

                int16_t qblk[8][8]; zz_to_block(zz,qblk);

                const int (*qt)[8]=(comp==0)?qtY:qtC;
                double F[8][8];
                for(int u=0;u<8;u++) for(int v=0;v<8;v++) F[u][v]=(double)qblk[u][v]*(double)qt[u][v];

                double blk[8][8];
                idct8x8(F,blk);

                int16_t *plane=(comp==0)?ycc.Y:(comp==1)?ycc.Cb:ycc.Cr;
                put_block_i16(plane,(int)W,(int)H,(int)bx,(int)by,blk);
            }
        }
    }

    fclose(fp);
    ImageRGB rgb=ycc_to_rgb_rct(&ycc);
    bmp_write_24(outbmp,&rgb);
    img_free(&rgb);
    ycc_free(&ycc);
}

// ---------- Huffman decode (Mode 3) ----------
typedef struct HNode {
    int leaf;
    uint32_t sym;
    struct HNode *l,*r;
} HNode;

static HNode* hn_new(void){ HNode *n=(HNode*)calloc(1,sizeof(HNode)); if(!n) die("OOM"); return n; }
static void hn_free(HNode *n){ if(!n) return; hn_free(n->l); hn_free(n->r); free(n); }

static void insert_code(HNode *root, const char *bits, uint32_t sym){
    HNode *cur=root;
    for(const char *p=bits; *p; p++){
        if(*p=='\n'||*p=='\r') break;
        if(*p=='0'){
            if(!cur->l) cur->l=hn_new();
            cur=cur->l;
        } else if(*p=='1'){
            if(!cur->r) cur->r=hn_new();
            cur=cur->r;
        }
    }
    cur->leaf=1;
    cur->sym=sym;
}

static HNode* build_tree_from_codebook(const char *codebook){
    FILE *fp=fopen(codebook,"r"); if(!fp) die("open codebook.txt");
    char line[4096];
    HNode *root=hn_new();
    while(fgets(line,sizeof(line),fp)){
        if(line[0]=='#' || line[0]=='\n') continue;
        unsigned skip; int val; unsigned cnt;
        char bits[2048]={0};
        if(sscanf(line,"%u %d %u %2047s",&skip,&val,&cnt,bits)==4){
            uint32_t sym=((uint32_t)skip<<16) | (uint16_t)(int16_t)val;
            insert_code(root,bits,sym);
        }
    }
    fclose(fp);
    return root;
}

static int next_symbol_from_ascii(FILE *fp, HNode *root, uint32_t *sym){
    char line[4096];
    if(!fgets(line,sizeof(line),fp)) return 0;
    HNode *cur=root;
    for(char *p=line; *p; p++){
        if(*p=='0') cur=cur->l;
        else if(*p=='1') cur=cur->r;
        else if(*p=='\n'||*p=='\r') break;
        if(!cur) die("bad huffman codeword");
    }
    if(!cur->leaf) die("non-leaf codeword");
    *sym=cur->sym;
    return 1;
}

typedef struct {
    FILE *fp;
    uint8_t cur;
    int bits_left;
    uint64_t nbits_total;
    uint64_t nbits_used;
} BitReader;

static void br_init(BitReader *br, FILE *fp, uint64_t nbits){
    br->fp=fp; br->cur=0; br->bits_left=0; br->nbits_total=nbits; br->nbits_used=0;
}
static int br_getbit(BitReader *br, int *bit){
    if(br->nbits_used >= br->nbits_total) return 0;
    if(br->bits_left==0){
        int c=fgetc(br->fp);
        if(c==EOF) die("unexpected EOF in bitstream");
        br->cur=(uint8_t)c;
        br->bits_left=8;
    }
    int b=(br->cur & 0x80)?1:0;
    br->cur <<= 1;
    br->bits_left--;
    br->nbits_used++;
    *bit=b;
    return 1;
}

static int next_symbol_from_bits(BitReader *br, HNode *root, uint32_t *sym){
    HNode *cur=root;
    while(cur && !cur->leaf){
        int b;
        if(!br_getbit(br,&b)) return 0;
        cur = (b==0)?cur->l:cur->r;
    }
    if(!cur) die("invalid huffman path");
    *sym=cur->sym;
    return 1;
}

static void mode3_ascii(const char *outbmp, const char *codebook, const char *hufftxt){
    HNode *root=build_tree_from_codebook(codebook);
    FILE *fp=fopen(hufftxt,"r"); if(!fp) die("open huffman_code.txt");

    int W,H;
    if(fscanf(fp,"%d %d",&W,&H)!=2) die("huffman txt size row");
    char tmp[256]; fgets(tmp,sizeof(tmp),fp);

    static const int qtY[8][8]={
     {16,11,10,16,24,40,51,61},{12,12,14,19,26,58,60,55},{14,13,16,24,40,57,69,56},{14,17,22,29,51,87,80,62},
     {18,22,37,56,68,109,103,77},{24,35,55,64,81,104,113,92},{49,64,78,87,103,121,120,101},{72,92,95,98,112,100,103,99}
    };
    static const int qtC[8][8]={
     {17,18,24,47,99,99,99,99},{18,21,26,66,99,99,99,99},{24,26,56,99,99,99,99,99},{47,66,99,99,99,99,99,99},
     {99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99}
    };

    int nbx=(W+7)/8, nby=(H+7)/8;

    ImageYCCi ycc; ycc.w=W; ycc.h=H;
    size_t n=(size_t)W*H;
    ycc.Y =(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cb=(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cr=(int16_t*)calloc(n,sizeof(int16_t));
    if(!ycc.Y||!ycc.Cb||!ycc.Cr) die("OOM");

    int16_t prevDC[3]={0,0,0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                int16_t zz[64]={0};
                int k=0;
                while(k<64){
                    uint32_t sym;
                    if(!next_symbol_from_ascii(fp,root,&sym)) die("EOF huffman symbols");
                    uint8_t skip=(uint8_t)(sym>>16);
                    int16_t val =(int16_t)(uint16_t)(sym & 0xFFFFu);
                    k += (int)skip;
                    if(k>=64) break;
                    zz[k]=val;
                    k++;
                }

                inv_dpcm_dc(zz,&prevDC[comp]);

                int16_t qblk[8][8]; zz_to_block(zz,qblk);
                const int (*qt)[8]=(comp==0)?qtY:qtC;

                double F[8][8];
                for(int u=0;u<8;u++) for(int v=0;v<8;v++) F[u][v]=(double)qblk[u][v]*(double)qt[u][v];

                double blk[8][8];
                idct8x8(F,blk);

                int16_t *plane=(comp==0)?ycc.Y:(comp==1)?ycc.Cb:ycc.Cr;
                put_block_i16(plane,W,H,bx,by,blk);
            }
        }
    }

    fclose(fp);
    hn_free(root);

    ImageRGB rgb=ycc_to_rgb_rct(&ycc);
    bmp_write_24(outbmp,&rgb);
    img_free(&rgb);
    ycc_free(&ycc);
}

static void mode3_binary(const char *outbmp, const char *codebook, const char *huffbin){
    HNode *root=build_tree_from_codebook(codebook);
    FILE *fp=fopen(huffbin,"rb"); if(!fp) die("open huffman_code.bin");

    uint32_t W,H;
    uint64_t nbits;
    if(fread(&W,4,1,fp)!=1) die("read W");
    if(fread(&H,4,1,fp)!=1) die("read H");
    if(fread(&nbits,8,1,fp)!=1) die("read nbits");

    static const int qtY[8][8]={
     {16,11,10,16,24,40,51,61},{12,12,14,19,26,58,60,55},{14,13,16,24,40,57,69,56},{14,17,22,29,51,87,80,62},
     {18,22,37,56,68,109,103,77},{24,35,55,64,81,104,113,92},{49,64,78,87,103,121,120,101},{72,92,95,98,112,100,103,99}
    };
    static const int qtC[8][8]={
     {17,18,24,47,99,99,99,99},{18,21,26,66,99,99,99,99},{24,26,56,99,99,99,99,99},{47,66,99,99,99,99,99,99},
     {99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99},{99,99,99,99,99,99,99,99}
    };

    int nbx=((int)W+7)/8, nby=((int)H+7)/8;

    ImageYCCi ycc; ycc.w=(int)W; ycc.h=(int)H;
    size_t n=(size_t)W*H;
    ycc.Y =(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cb=(int16_t*)calloc(n,sizeof(int16_t));
    ycc.Cr=(int16_t*)calloc(n,sizeof(int16_t));
    if(!ycc.Y||!ycc.Cb||!ycc.Cr) die("OOM");

    BitReader br; br_init(&br,fp,nbits);
    int16_t prevDC[3]={0,0,0};

    for(int by=0;by<nby;by++){
        for(int bx=0;bx<nbx;bx++){
            for(int comp=0;comp<3;comp++){
                int16_t zz[64]={0};
                int k=0;
                while(k<64){
                    uint32_t sym;
                    if(!next_symbol_from_bits(&br,root,&sym)) die("EOF bitstream too early");
                    uint8_t skip=(uint8_t)(sym>>16);
                    int16_t val =(int16_t)(uint16_t)(sym & 0xFFFFu);
                    k += (int)skip;
                    if(k>=64) break;
                    zz[k]=val;
                    k++;
                }

                inv_dpcm_dc(zz,&prevDC[comp]);

                int16_t qblk[8][8]; zz_to_block(zz,qblk);
                const int (*qt)[8]=(comp==0)?qtY:qtC;

                double F[8][8];
                for(int u=0;u<8;u++) for(int v=0;v<8;v++) F[u][v]=(double)qblk[u][v]*(double)qt[u][v];

                double blk[8][8];
                idct8x8(F,blk);

                int16_t *plane=(comp==0)?ycc.Y:(comp==1)?ycc.Cb:ycc.Cr;
                put_block_i16(plane,(int)W,(int)H,bx,by,blk);
            }
        }
    }

    fclose(fp);
    hn_free(root);

    ImageRGB rgb=ycc_to_rgb_rct(&ycc);
    bmp_write_24(outbmp,&rgb);
    img_free(&rgb);
    ycc_free(&ycc);
}

// ---------- Mode 1 wrapper: decide 1(a) vs 1(b) by argc ----------
static void mode1_dispatch(int argc, char **argv){
    // 1(a): decoder 1 QRes.bmp orig.bmp QtY QtCb QtCr dim qFY qFCb qFCr   (argc == 11)
    // 1(b): decoder 1 Res.bmp QtY QtCb QtCr dim qFY qFCb qFCr eFY eFCb eFCr (argc == 12)
    if(argc!=11 && argc!=12) die("decoder 1: wrong args (see spec)");

    const char *outbmp=argv[2];
    int use_error = (argc==12);

    const char *origbmp = use_error ? NULL : argv[3];

    const char *QtYp = use_error ? argv[3] : argv[4];
    const char *QtCbp= use_error ? argv[4] : argv[5];
    const char *QtCrp= use_error ? argv[5] : argv[6];
    const char *dimp = use_error ? argv[6] : argv[7];

    const char *qFYp = use_error ? argv[7] : argv[8];
    const char *qFCbp= use_error ? argv[8] : argv[9];
    const char *qFCrp= use_error ? argv[9] : argv[10];

    const char *eFYp = use_error ? argv[10] : NULL;
    const char *eFCbp= use_error ? argv[11] : NULL;
    const char *eFCrp= use_error ? argv[12] : NULL; // won't be used; argc==12 means argv[11] exists; but this line is safe only if argc>=13 (not)
    // Fix indexing for argc==12:
    if(use_error){
        eFYp = argv[10];
        eFCbp= argv[11];
        eFCrp= argv[12-1]; // argv[11] is eFCb, argv[11]?? -> Actually argc==12 => argv[11] exists, argv[12] doesn't.
        // We'll map correctly below with explicit.
    }

    // parse dim
    int W,H;
    FILE *fd=fopen(dimp,"r"); if(!fd) die("open dim.txt");
    if(fscanf(fd,"%d %d",&W,&H)!=2) die("dim parse");
    fclose(fd);

    int qtY[8][8], qtCb[8][8], qtCr[8][8];
    read_qt_txt(QtYp,qtY);
    read_qt_txt(QtCbp,qtCb);
    read_qt_txt(QtCrp,qtCr); (void)qtCr;

    FILE *fqY=fopen(qFYp,"rb"); if(!fqY) die("open qF_Y.raw");
    FILE *fqCb=fopen(qFCbp,"rb"); if(!fqCb) die("open qF_Cb.raw");
    FILE *fqCr=fopen(qFCrp,"rb"); if(!fqCr) die("open qF_Cr.raw");

    FILE *feY=NULL,*feCb=NULL,*feCr=NULL;

    if(use_error){
        // argc==12 layout:
        // decoder 1 Res.bmp QtY QtCb QtCr dim qFY qFCb qFCr eFY eFCb eFCr
        // argv indexes: 0..11
        const char *eFY = argv[10];
        const char *eFCb= argv[11];
        const char *eFCr= argv[11]; // placeholder, will be fixed by requiring argc==13 ideally
        // We can't know eFCr without having it in argv; so we enforce correct argc for 1(b):
        die("decoder 1(b): spec needs 11 args after '1' (total argc should be 13). Please call exactly like spec.");
    }

    // For correctness, we support exactly the spec calls:
    // 1(a) total argc = 11
    // 1(b) total argc = 12 in your earlier text includes 11 inputs after '1' => actually total argc should be 13.
    // So we only run 1(a) here; 1(b) implemented below with correct argc in main.

    decode_mode1(outbmp,W,H,qtY,qtCb,fqY,fqCb,fqCr,NULL,NULL,NULL,0);
    fclose(fqY); fclose(fqCb); fclose(fqCr);

    ImageRGB orig=bmp_read_24(origbmp);
    ImageRGB rec =bmp_read_24(outbmp);
    print_pixel_sqnr_rgb(&orig,&rec);
    img_free(&orig); img_free(&rec);
}

// ---------- main ----------
int main(int argc, char **argv){
    if(argc<2) die("usage: decoder <mode> ...");

    int mode=atoi(argv[1]);

    if(mode==0){
        if(argc!=7) die("decoder 0 out.bmp R.txt G.txt B.txt dim.txt");
        mode0(argv[2],argv[3],argv[4],argv[5],argv[6]);
        return 0;
    }

    if(mode==1){
        // Support spec EXACTLY:
        // 1(a): decoder 1 QRes.bmp Kimberly.bmp QtY QtCb QtCr dim qFY qFCb qFCr  => argc = 11
        // 1(b): decoder 1 Res.bmp QtY QtCb QtCr dim qFY qFCb qFCr eFY eFCb eFCr => argc = 12? actually: out + 10 inputs => total argc=13
        if(argc==11){
            mode1_dispatch(argc,argv);
            return 0;
        } else if(argc==12 || argc==13){
            // implement 1(b) with correct indexing if argc==13
            if(argc!=13) die("decoder 1(b) needs: decoder 1 Res.bmp QtY QtCb QtCr dim qFY qFCb qFCr eFY eFCb eFCr (total argc=13)");
            const char *outbmp=argv[2];
            const char *QtYp=argv[3], *QtCbp=argv[4], *QtCrp=argv[5], *dimp=argv[6];
            const char *qFYp=argv[7], *qFCbp=argv[8], *qFCrp=argv[9];
            const char *eFYp=argv[10], *eFCbp=argv[11], *eFCrp=argv[12];

            int W,H;
            FILE *fd=fopen(dimp,"r"); if(!fd) die("open dim.txt");
            if(fscanf(fd,"%d %d",&W,&H)!=2) die("dim parse");
            fclose(fd);

            int qtY[8][8], qtCb[8][8], qtCr[8][8];
            read_qt_txt(QtYp,qtY);
            read_qt_txt(QtCbp,qtCb);
            read_qt_txt(QtCrp,qtCr); (void)qtCr;

            FILE *fqY=fopen(qFYp,"rb"); if(!fqY) die("open qF_Y.raw");
            FILE *fqCb=fopen(qFCbp,"rb"); if(!fqCb) die("open qF_Cb.raw");
            FILE *fqCr=fopen(qFCrp,"rb"); if(!fqCr) die("open qF_Cr.raw");

            FILE *feY=fopen(eFYp,"rb"); if(!feY) die("open eF_Y.raw");
            FILE *feCb=fopen(eFCbp,"rb"); if(!feCb) die("open eF_Cb.raw");
            FILE *feCr=fopen(eFCrp,"rb"); if(!feCr) die("open eF_Cr.raw");

            decode_mode1(outbmp,W,H,qtY,qtCb,fqY,fqCb,fqCr,feY,feCb,feCr,1);

            fclose(fqY); fclose(fqCb); fclose(fqCr);
            fclose(feY); fclose(feCb); fclose(feCr);
            return 0;
        } else {
            die("decoder 1: wrong args");
        }
    }

    if(mode==2){
        if(argc!=5) die("decoder 2 out.bmp ascii|binary rle_code.(txt|bin)");
        if(strcmp(argv[3],"ascii")==0) mode2_ascii(argv[2],argv[4]);
        else if(strcmp(argv[3],"binary")==0) mode2_binary(argv[2],argv[4]);
        else die("mode2 fmt must be ascii or binary");
        return 0;
    }

    if(mode==3){
        if(argc!=6) die("decoder 3 out.bmp ascii|binary codebook.txt huffman_code.(txt|bin)");
        if(strcmp(argv[3],"ascii")==0) mode3_ascii(argv[2],argv[4],argv[5]);
        else if(strcmp(argv[3],"binary")==0) mode3_binary(argv[2],argv[4],argv[5]);
        else die("mode3 fmt must be ascii or binary");
        return 0;
    }

    die("mode must be 0..3");
    return 0;
}
