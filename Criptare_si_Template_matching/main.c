#include <stdio.h>
#include <stdlib.h>

/// ////// FUNCTII PENTRU CRIPTARE /////////////

///xorshift32
int xorshift32(unsigned int *r)
{

    *r=*r^*r<<13;
    *r=*r^*r>>17;
    *r=*r^*r<<5;
    return *r;

}


///incarcarea pixelilor imaginii BMP
int IncarcareBMP(char *s,int *n)
{
    unsigned char *L;
    FILE* fin=fopen(s,"rb");

    unsigned int width_img, height_img;
    fseek(fin, 18, SEEK_SET);
    fread(&width_img, sizeof(unsigned int), 1, fin);
    fread(&height_img, sizeof(unsigned int), 1, fin);
    unsigned int padding=0;

    if(width_img % 4)
        padding=4-(3*width_img)%4;

    L=(unsigned char*)malloc((width_img*3+padding)*height_img);

    fseek(fin,54,SEEK_SET);
    int k=-1;
    unsigned char c;
    while(fread(&c,1,1,fin)==1)
    {
        L[++k]=c;

    }
    *n=k+1;
    fclose(fin);
    return L;


}


///incarcarea headerului imaginii BMP
int IncarcareHeader(char *s)
{

    unsigned char *header;
    FILE*fin=fopen(s,"rb");
    header=(unsigned char*)malloc(54);
    int i;
    for(i=0; i<54; i++)
        fread(&header[i],1,1,fin);
    fclose(fin);
    return header;


}


///dezalocare vectori folositi in criptare
void dezalocareC(unsigned char *L,unsigned char *header,unsigned int *R,int *p,unsigned char *L1,unsigned char *C)
{
    free(L);
    free(header);
    free(R);
    free(p);
    free(L1);
    free(C);
}

///dezalocare vectori folositi in decriptare
void dezalocareD(unsigned char *C,unsigned char *header,unsigned int *R,int *p,int *p1,unsigned char *C1,unsigned char *D)
{
    free(C);
    free(header);
    free(R);
    free(p);
    free(p1);
    free(C1);
    free(D);
}

///incarcarea in memoria externa a noii imagini BMP
void DescarcareBMP(char *s,unsigned char *L,int n,unsigned char *header)
{
    FILE*fout=fopen(s,"wb");
    int i;

    for(i=0; i<54; i++)
        fwrite(&header[i],1,1,fout);

    for(i=0; i<n; i++)
        fwrite(&L[i],1,1,fout);


    fclose(fout);
}


///criptarea propriu-zisa
void criptare(char *s1,char *s2,char *s3)
{
    ///fisierele necesare criptarii
    FILE*fin=fopen(s1,"rb");
    FILE*fout=fopen(s2,"wb");
    FILE*fkey=fopen(s3,"r");

    ///vectorii de care va fi nevoie pentru a crea imaginea criptata
    unsigned char *L;
    unsigned char *header;
    int n;

    L=IncarcareBMP(s1,&n);
    header=IncarcareHeader(s1);

    ///dimensiunile imaginii initiale
    unsigned int width_img, height_img;
    fseek(fin, 18, SEEK_SET);
    fread(&width_img, sizeof(unsigned int), 1, fin);
    fread(&height_img, sizeof(unsigned int), 1, fin);

    ///valorile din fisierul cu cheia secreta
    unsigned int R0,SV;
    fscanf(fkey,"%u",&R0);
    fscanf(fkey,"%u",&SV);


    ///a) sirul de numere pseudo-aleatoare
    unsigned int *R;
    R=(unsigned int*)malloc(2*width_img*height_img*sizeof(unsigned int));

    int i;
    R[0]=R0;
    for(i=1; i<=2*width_img*height_img-1; i++)
        R[i]=xorshift32(&R[i-1]);

    ///b)permutarea aleatoare
    int r,*p,aux;

    p=(int*)malloc(width_img*height_img*sizeof(int));

    for(i=0; i<width_img*height_img; i++)
        p[i]=i;

    for(i=width_img*height_img-1; i>=1; i--)
    {
        r=R[i]%i;
        aux=p[r];
        p[r]=p[i];
        p[i]=aux;
    }

    ///permutarea pixelilor (grupe de cate 3 octeti)
    unsigned char *L1;

    unsigned int padding=0;
    if(width_img % 4)
        padding=4-(3*width_img)%4;

    L1=(unsigned char*)malloc((width_img*3+padding)*height_img);

    for(i=0; i<width_img*height_img; i++)
    {
        L1[p[i]*3]=L[i*3];
        L1[p[i]*3+1]=L[i*3+1];
        L1[p[i]*3+2]=L[i*3+2];
    }

///schimbarea valorilor pixelilor
    unsigned char *C;
    C=(unsigned char*)malloc((width_img*3+padding)*height_img);

    unsigned char pix[3],nr1[4],nr2[4];

    for(i=0; i<3*width_img*height_img; i+=3)
    {
        if(i==0)
        {
            pix[0]=L1[0];
            pix[1]=L1[1];
            pix[2]=L1[2];
            unsigned int save;
            save=SV;
            save=save>>24;
            nr1[3]=save;
            save=SV;
            save=save<<8;
            save=save>>24;
            nr1[2]=save;
            save=SV;
            save=save<<16;
            save=save>>24;
            nr1[1]=save;
            save=SV;
            save=save<<24;
            save=save>>24;
            nr1[0]=save;

            save=R[width_img*height_img];
            save=save>>24;
            nr2[3]=save;
            save=R[width_img*height_img];
            save=save<<8;
            save=save>>24;
            nr2[2]=save;
            save=R[width_img*height_img];
            save=save<<16;
            save=save>>24;
            nr2[1]=save;
            save=R[width_img*height_img];
            save=save<<24;
            save=save>>24;
            nr2[0]=R[width_img*height_img];

            C[i]=(pix[0]^nr1[2])^nr2[2];
            C[i+1]=(pix[1]^nr1[1])^nr2[1];
            C[i+2]=(pix[2]^nr1[0])^nr2[0];

        }

        else
        {
            pix[0]=L1[i];
            pix[1]=L1[i+1];
            pix[2]=L1[i+2];

            unsigned int save=R[width_img*height_img+(i/3)];
            save=save>>24;
            nr1[3]=save;
            save=R[width_img*height_img+(i/3)];
            save=save<<8;
            save=save>>24;
            nr1[2]=save;
            save=R[width_img*height_img+(i/3)];
            save=save<<16;
            save=save>>24;
            nr1[1]=save;
            save=R[width_img*height_img+(i/3)];
            save=save<<24;
            save=save>>24;
            nr1[0]=save;

            C[i]=(pix[0]^C[i-3])^nr1[2];
            C[i+1]=(pix[1]^C[i-2])^nr1[1];
            C[i+2]=(pix[2]^C[i-1])^nr1[0];

        }



    }
///crearea imaginii criptate si dezalocarea memoriei
    DescarcareBMP(s2,C,n,header);
    dezalocareC(L,header,R,p,L1,C);

    fclose(fin);
    fclose(fout);
    fclose(fkey);

}


///decriptrea
void decriptare(char *s2,char *s3,char *s4)
{
    ///fisierele necesare
    FILE*fin=fopen(s2,"rb");
    FILE*fout=fopen(s3,"wb");
    FILE*fkey=fopen(s4,"r");

    ///vectorul ce contile pixelii imaginii criptate
    unsigned char *C;
    unsigned char *header;
    int n;

    C=IncarcareBMP(s2,&n);
    header=IncarcareHeader(s2);

    ///dimensiunile imaginii criptate
    unsigned int width_img, height_img;
    fseek(fin, 18, SEEK_SET);
    fread(&width_img, sizeof(unsigned int), 1, fin);
    fread(&height_img, sizeof(unsigned int), 1, fin);

    unsigned int padding=0;
    if(width_img % 4)
        padding=4-(3*width_img)%4;

    ///valorile din fisierul cu cheia secreta
    unsigned int R0,SV;
    fscanf(fkey,"%u",&R0);
    fscanf(fkey,"%u",&SV);


    ///a) sirul de numere pseudo-aleatoare
    unsigned int *R;
    R=(unsigned int*)malloc(2*width_img*height_img*sizeof(unsigned int));

    int i;
    R[0]=R0;
    for(i=1; i<=2*width_img*height_img-1; i++)
        R[i]=xorshift32(&R[i-1]);

    ///b)permutarea aleatoare
    int r,*p,aux;

    p=(int*)malloc(width_img*height_img*sizeof(int));

    for(i=0; i<width_img*height_img; i++)
        p[i]=i;

    for(i=width_img*height_img-1; i>=1; i--)
    {
        r=R[i]%i;
        aux=p[r];
        p[r]=p[i];
        p[i]=aux;
    }

    int *p1;///inversa permutarii aleatoare p
    p1=(int*)malloc(width_img*height_img*sizeof(int));

    for(i=0; i<width_img*height_img; i++)
        p1[p[i]]=i;

    ///schimbarea valorilor octetilor din C (vectorul in care sunt pixelii criptati)
    unsigned char *C1;
    C1=(unsigned char*)malloc((width_img*3+padding)*height_img);

    unsigned char pix[3],nr1[4],nr2[4];

    for(i=0; i<3*width_img*height_img; i+=3)
    {
        if(i==0)
        {
            pix[0]=C[0];
            pix[1]=C[1];
            pix[2]=C[2];
            unsigned int save;
            save=SV;
            save=save>>24;
            nr1[3]=save;
            save=SV;
            save=save<<8;
            save=save>>24;
            nr1[2]=save;
            save=SV;
            save=save<<16;
            save=save>>24;
            nr1[1]=save;
            save=SV;
            save=save<<24;
            save=save>>24;
            nr1[0]=save;

            save=R[width_img*height_img];
            save=save>>24;
            nr2[3]=save;
            save=R[width_img*height_img];
            save=save<<8;
            save=save>>24;
            nr2[2]=save;
            save=R[width_img*height_img];
            save=save<<16;
            save=save>>24;
            nr2[1]=save;
            save=R[width_img*height_img];
            save=save<<24;
            save=save>>24;
            nr2[0]=R[width_img*height_img];

            C1[i]=(pix[0]^nr1[2])^nr2[2];
            C1[i+1]=(pix[1]^nr1[1])^nr2[1];
            C1[i+2]=(pix[2]^nr1[0])^nr2[0];

        }

        else
        {
            pix[0]=C[i];
            pix[1]=C[i+1];
            pix[2]=C[i+2];

            unsigned int save=R[width_img*height_img+(i/3)];
            save=save>>24;
            nr1[3]=save;
            save=R[width_img*height_img+(i/3)];
            save=save<<8;
            save=save>>24;
            nr1[2]=save;
            save=R[width_img*height_img+(i/3)];
            save=save<<16;
            save=save>>24;
            nr1[1]=save;
            save=R[width_img*height_img+(i/3)];
            save=save<<24;
            save=save>>24;
            nr1[0]=save;

            C1[i]=(pix[0]^C[i-3])^nr1[2];
            C1[i+1]=(pix[1]^C[i-2])^nr1[1];
            C1[i+2]=(pix[2]^C[i-1])^nr1[0];

        }

    }

///permutarea pixelilor si obtinerea imaginii decriptate
    unsigned char *D;
    D=(unsigned char*)malloc((width_img*3+padding)*height_img);

    for(i=0; i<width_img*height_img; i++)
    {
        D[p1[i]*3]=C1[i*3];
        D[p1[i]*3+1]=C1[i*3+1];
        D[p1[i]*3+2]=C1[i*3+2];
    }

///salvarea in memoria externa a imaginii decriptate(initiale) si dezalocarea memoriei
    DescarcareBMP(s3,D,n,header);
    dezalocareD(C,header,R,p,p1,C1,D);





}

///functia care afiseaza pe ecran valorile testului chi-patrat pentru o imagine introdusa
void testul_chi_patrat(char *s)
{
    FILE*fin=fopen(s,"rb");

    unsigned char *V;
    int n;
    V=IncarcareBMP(s,&n);


    unsigned int width_img, height_img;
    fseek(fin, 18, SEEK_SET);
    fread(&width_img, sizeof(unsigned int), 1, fin);
    fread(&height_img, sizeof(unsigned int), 1, fin);

    double f;
    f=(width_img*height_img)/256;

    int i,j;
    double XR=0,XG=0,XB=0;

    for(i=0; i<=255; i++)
    {
        int frecv=0;
        for(j=0; j<n; j+=3)
            if(V[j]==i)
                frecv++;
        XR+=((frecv-f)*(frecv-f))/f;

    }

    for(i=0; i<=255; i++)
    {
        int frecv=0;
        for(j=1; j<n; j+=3)
            if(V[j]==i)
                frecv++;
        XG+=((frecv-f)*(frecv-f))/f;

    }

    for(i=0; i<=255; i++)
    {
        int frecv=0;
        for(j=2; j<n; j+=3)
            if(V[j]==i)
                frecv++;
        XB+=((frecv-f)*(frecv-f))/f;

    }

    printf("%lf\n%lf\n%lf\n\n",XR,XG,XB);

    free(V);
    fclose(fin);




}


/// ////// FUNCTII PENTRU TEMPLATE MATCHING /////////////

///structura cu ajutorul careia vor fi salvate culorile pentru detectii
struct culoare
{
    unsigned char r;
    unsigned char g;
    unsigned char b;
};


///functia care transforma o imagine in imagine grayscale
void grayscale(char* s1,char* s2)
{
    FILE *fin, *fout;
    unsigned int width_img, height_img;
    unsigned char pixel[3], header[54], aux;


    fin = fopen(s1, "rb");
    fout = fopen(s2, "wb+");

    fseek(fin, 18, SEEK_SET);
    fread(&width_img, sizeof(unsigned int), 1, fin);
    fread(&height_img, sizeof(unsigned int), 1, fin);

    fseek(fin,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,fin)==1)
    {
        fwrite(&c,1,1,fout);
        fflush(fout);
    }
    fclose(fin);


    int padding=0;
    if(width_img % 4 != 0)
        padding = 4 - (3 * width_img) % 4;

    fseek(fout, 54, SEEK_SET);
    int i,j;
    for(i = 0; i < height_img; i++)
    {
        for(j = 0; j < width_img; j++)
        {
            fread(pixel, 3, 1, fout);
            aux = 0.299*pixel[2] + 0.587*pixel[1] + 0.114*pixel[0];
            pixel[0] = pixel[1] = pixel[2] = aux;
            fseek(fout, -3, SEEK_CUR);
            fwrite(pixel, 3, 1, fout);
            fflush(fout);
        }
        fseek(fout,padding,SEEK_CUR);
    }
    fclose(fout);
}

///functia care incarca o imagine intr-o matrice de care va fi nevoie in functa template_matching
int incarcare_imagine(char *s,unsigned char **m,unsigned int *w,unsigned int *h)
{
    FILE*fin=fopen(s,"rb");

    fseek(fin,18,SEEK_SET);
    fread(&(*w),sizeof(unsigned int),1,fin);
    fread(&(*h),sizeof(unsigned int),1,fin);

    int padding=0;
    if(*w % 4 != 0)
        padding = 4 - (3 * *w) % 4;

    m=(unsigned char**)malloc((*h*sizeof(unsigned char*)));
    int i,j;

    for(i=0; i<*h; i++)
        m[i]=(unsigned char*)malloc((*w*3+padding)*sizeof(unsigned char));///fara padding


    fseek(fin,54,SEEK_SET);
    for(i=0; i<*h; i++)
        for(j=0; j<(*w)*3+padding; j++)
            fread(&m[i][j],1,1,fin);

    fclose(fin);
    return m;

}

///functia care incarca header-ul imaginii, care va fi folosit in constuirea noilor imagini
int incarcare_header(char *s)
{

    unsigned char *header;
    FILE*fin=fopen(s,"rb");
    header=(unsigned char*)malloc(54);
    int i;
    for(i=0; i<54; i++)
        fread(&header[i],1,1,fin);
    fclose(fin);
    return header;


}

///functia care construieste in memoria externa o imagine cu ajutorul unei matrice si al header-ului
void descarcare_imagine(char *s,unsigned char **m,unsigned int w,unsigned int h,unsigned char *header)
{

    FILE*fout=fopen(s,"wb");

    int i,j;

    for(i=0; i<54; i++)
        fwrite(&header[i],1,1,fout);

    for(i=0; i<h; i++)
        for(j=0; j<w*3; j++)
            fwrite(&m[i][j],1,1,fout);

    fclose(fout);

}



///funcita propriu-zisa care gliseaza sablonul pe imaginea initiala si determina corelatiile
void template_matching(unsigned char **I,unsigned char **T,unsigned int wi,unsigned int hi,double ps,char *sS,struct culoare clr[],int nr_sablon,double **corr)
{

    unsigned char **S;
    unsigned int ws;
    unsigned int hs;

///incarcarea sablonului
    S=incarcare_imagine(sS,S,&ws,&hs);

    int i,j;

///initializarea matricei de corelatii cu 0
    for(i=0; i<hi-hs; i++)
        for(j=0; j<wi-ws; j++)
            corr[i][j]=0;



    unsigned int wc=wi-ws,hc=hi-hs;

///variabile utilizare la determinarea corelatiilor
    int ns=ws*hs;
    double fi,si;
    int nfi,nsi;
    double Ds,Dfi;

///calcularea corelatiilor si salvarea acestora in matrice
    int k,l;
    for(i=0; i<hc; i++)
        for(j=0; j<3*wc; j+=3)
        {
            fi=0;
            nfi=0;
            for(k=i; k<hs+i; k++)
                for(l=j; l<3*ws+j; l++)
                {
                    fi+=I[k][l];
                    nfi++;
                }
            fi/=nfi;

            si=0;
            nsi=0;
            for(k=0; k<hs; k++)
                for(l=0; l<3*ws; l++)
                {
                    si+=S[k][l];
                    nsi++;
                }
            si/=nsi;



            Dfi=0;
            for(k=i; k<hs+i; k++)
                for(l=j; l<3*ws+j; l+=3)
                    Dfi+=(I[k][l]-fi)*(I[k][l]-fi);
            Dfi/=ns-1;
            Dfi=sqrt(Dfi);

            Ds=0;
            for(k=0; k<hs; k++)
                for(l=0; l<3*ws; l++)
                    Ds+=(S[k][l]-si)*(S[k][l]-si);
            Ds/=ns-1;
            Ds=sqrt(Ds);

            for(k=i; k<hs+i; k++)
                for(l=j; l<3*ws+j; l+=3)
                    corr[i][j]+=(1/(Dfi*Ds))*(I[k][l]-fi)*(S[k-i][l-j]-si);

            corr[i][j]/=(double)ns;

        }

///conturarea cifrelor recunoscute, a caror corelatie este mai mare decat pragul ps
    for(i=0; i<hc; i++)
        for(j=0; j<3*wc; j+=3)
            if(corr[i][j]>=ps)
                desenare_contur(T,i,j,clr,nr_sablon,hs,ws);

///dezalocare memorie
    for(i=0; i<hs; i++)
        free(S[i]);
    free(S);



}


///functia care deseneaza pe poza initiala conturul cifrei recunoscute incepand de la pozitia (i,j) din matricea imaginii
void desenare_contur(unsigned char **I,int i,int j,struct culoare clr[],int nrc,int hs,int ws)
{
    int k,l;

    ///linia din stanga
    for(k=i; k<hs+i; k++)
    {
        I[k][j]=clr[nrc].b;
        I[k][j+1]=clr[nrc].g;
        I[k][j+2]=clr[nrc].r;

    }

    ///linia de sus
    for(l=j; l<3*ws+j; l+=3)
    {
        I[i][l]=clr[nrc].b;
        I[i][l+1]=clr[nrc].g;
        I[i][l+2]=clr[nrc].r;
    }

    ///linia din dreapta
    for(k=i; k<hs+i; k++)
    {
        I[k][3*ws+j]=clr[nrc].b;
        I[k][3*ws+j+1]=clr[nrc].g;
        I[k][3*ws+j+2]=clr[nrc].r;

    }

    ///linia de jos
    for(l=j; l<3*ws+j; l+=3)
    {
        I[hs+i-1][l]=clr[nrc].b;
        I[hs+i-1][l+1]=clr[nrc].g;
        I[hs+i-1][l+2]=clr[nrc].r;
    }
}

///functia folosita in qsort pentru a sorta descrescator vectorul de corelatii
int sort(const void *a,const void *b)
{


    double i1 = *(double*) b;
    double i2 = *(double*) a;
    if (i1 < i2)
        return -1;
    else if (i1 == i2)
        return 0;
    else
        return 1;


}

///aceasta functie nu este testata sau aplicata..
void eliminare_nonmaxime(double *D,double **corr0,double **corr1,double **corr2,double **corr3,double **corr4,double **corr5,double **corr6,double **corr7,double **corr8,double **corr9,unsigned int w,unsigned int h,int nrd)
{
    int i,j,k,l;
    h-=15;
    w-=11;

    for(i=0; i<nrd-1; i++)
    {
        int i1=-1,j1=-1;

        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr0[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }

        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr1[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }

        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr2[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr3[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr4[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr5[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr6[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr7[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr8[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }
        if(i1==-1&&j1==-1)
            for(k=0; k<h; k++)
                for(l=0; l<w; l++)
                    if(corr9[k][l]==D[i])
                    {
                        i1=k;
                        j1=l;
                    }


        for(j=i+1; j<nrd; j++)
        {
            int i2=-1,j2=-1;


            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr0[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }

            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr1[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }

            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr2[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr3[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr4[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr5[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr6[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr7[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i1==-1&&j1==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr8[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }
            if(i2==-1&&j2==-1)
                for(k=0; k<h; k++)
                    for(l=0; l<w; l++)
                        if(corr9[k][l]==D[j])
                        {
                            i2=k;
                            j2=l;
                        }


            int mi,Mi,mj,Mj;
            if(i1>i2)
            {
                Mi=i1;
                mi=i2;
            }
            else
            {
                Mi=i2;
                mi=i1;
            }
            if(j1>j2)
            {
                Mj=j1;
                mj=j2;
            }
            else
            {
                Mj=j2;
                mj=j1;
            }

            double suprapunere;

            if(mi+15>Mi&&mj+11>Mj)

            {
                suprapunere=((Mi-mi)*(Mj-mj))/((15*11)+(15*11)-((Mi-mi)*(Mj-mj)));
                if  (suprapunere>0.2)
                    D[j]=-10;
            }



        }



    }

}





int main()
{

    char nume_secretkey[]="secret_key.txt";

    char nume_img1[] = "test.bmp";
    char nume_img1_criptat[]="p_criptat.bmp";
    char nume_img1_decriptat[]="p_decriptat.bmp";
    char nume_img2[] = "test_gs.bmp";
    char nume_img3[]="test_detectii.bmp";

    char nume_0[]="cifra0.bmp";
    char nume_gs0[]="cifra0_gs.bmp";
    char nume_1[]="cifra1.bmp";
    char nume_gs1[]="cifra1_gs.bmp";
    char nume_2[]="cifra2.bmp";
    char nume_gs2[]="cifra2_gs.bmp";
    char nume_3[]="cifra3.bmp";
    char nume_gs3[]="cifra3_gs.bmp";
    char nume_4[]="cifra4.bmp";
    char nume_gs4[]="cifra4_gs.bmp";
    char nume_5[]="cifra5.bmp";
    char nume_gs5[]="cifra5_gs.bmp";
    char nume_6[]="cifra6.bmp";
    char nume_gs6[]="cifra6_gs.bmp";
    char nume_7[]="cifra7.bmp";
    char nume_gs7[]="cifra7_gs.bmp";
    char nume_8[]="cifra8.bmp";
    char nume_gs8[]="cifra8_gs.bmp";
    char nume_9[]="cifra9.bmp";
    char nume_gs9[]="cifra9_gs.bmp";

/// PARTEA DE CRIPTARE

    ///criptarea unei imagini color BMP si salvarea acesteia in memoria externa

    criptare(nume_img1,nume_img1_criptat,nume_secretkey);


    ///decriptarea unei imagini color BMP si salvarea acesteia in memoria externa

    decriptare(nume_img1_criptat,nume_img1_decriptat,nume_secretkey);


    ///afisarea pe ecran a valorilor testului chi-patrat pentru imagina intiala si criptata, pe fiecare canal de culoare

    printf("Valorile testului chi-patrat pentru imaginea criptata:\n");
    testul_chi_patrat(nume_img1_criptat);
    printf("\nValorile testului chi-patrat pentru imaginea initiala/decriptata:\n");
    testul_chi_patrat(nume_img1_decriptat);


///PARTEA DE TEMPLATE MATCHING

///transformarea imaginii initiale si a sabloanelor in imagini grayscale
    grayscale(nume_img1,nume_img2);
    grayscale(nume_0,nume_gs0);
    grayscale(nume_1,nume_gs1);
    grayscale(nume_2,nume_gs2);
    grayscale(nume_3,nume_gs3);
    grayscale(nume_4,nume_gs4);
    grayscale(nume_5,nume_gs5);
    grayscale(nume_6,nume_gs6);
    grayscale(nume_7,nume_gs7);
    grayscale(nume_8,nume_gs8);
    grayscale(nume_9,nume_gs9);


    ///definirea culorilor

    struct culoare clr[10];

    ///rosu
    clr[0].r=255;
    clr[0].g=0;
    clr[0].b=0;

    ///galben
    clr[1].r=255;
    clr[1].g=255;
    clr[1].b=0;

    ///verde
    clr[2].r=0;
    clr[2].g=255;
    clr[2].b=0;

    ///cyan
    clr[3].r=0;
    clr[3].g=255;
    clr[3].b=255;

    ///magenta
    clr[4].r=255;
    clr[4].g=0;
    clr[4].b=255;

    ///albastru
    clr[5].r=0;
    clr[5].g=0;
    clr[5].b=255;

    ///argintiu
    clr[6].r=192;
    clr[6].g=192;
    clr[6].b=192;

    ///albastru
    clr[7].r=255;
    clr[7].g=140;
    clr[7].b=0;

    ///magenta
    clr[8].r=128;
    clr[8].g=0;
    clr[8].b=128;

    ///albastru
    clr[9].r=128;
    clr[9].g=0;
    clr[9].b=0;



    ///conturarea tuturor detectiilor

    unsigned char **I;///imaginea grayscale
    unsigned char **T;///imaginea initiala
    unsigned int W,H;///dimensiunile acestor imagini

    T=incarcare_imagine(nume_img1,T,&W,&H);
    I=incarcare_imagine(nume_img2,I,&W,&H);


    unsigned char *header;
    header=incarcare_header(nume_img2);


    ///declarrea matricelor de corelatie a tuturor celor 10 sabloane
    double **corr0,**corr1,**corr2,**corr3,**corr4,**corr5,**corr6,**corr7,**corr8,**corr9;

    int i,j;

///alocarea matricelor corelatiilor fiecarui sablon
    corr0=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr0[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr1=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr1[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr2=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr2[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr3=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr3[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr4=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr4[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr5=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr5[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr6=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr6[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr7=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr7[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr8=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr8[i]=(double*)malloc((W-11)*3*sizeof(double));

    corr9=(double**)malloc((H-15)*sizeof(double*));
    for(i=0; i<H-15; i++)
        corr9[i]=(double*)malloc((W-11)*3*sizeof(double));



///aplicarea operatiei de template matching asupra fiecarui sablon si salvarea corelatiilor acestuia in matricea corespunzatoare
    template_matching(I,T,W,H,0.35,nume_gs0,clr,0,corr0);
    template_matching(I,T,W,H,0.35,nume_gs1,clr,1,corr1);
    template_matching(I,T,W,H,0.35,nume_gs2,clr,2,corr2);
    template_matching(I,T,W,H,0.35,nume_gs3,clr,3,corr3);
    template_matching(I,T,W,H,0.35,nume_gs4,clr,4,corr4);
    template_matching(I,T,W,H,0.35,nume_gs5,clr,5,corr5);
    template_matching(I,T,W,H,0.35,nume_gs6,clr,6,corr6);
    template_matching(I,T,W,H,0.35,nume_gs7,clr,7,corr7);
    template_matching(I,T,W,H,0.35,nume_gs8,clr,8,corr8);
    template_matching(I,T,W,H,0.35,nume_gs9,clr,9,corr9);

    ///salvarea in memoria externa a imaginii pe care sunt conturate toate detectiile mai mari ca pragul introdus
    descarcare_imagine(nume_img3,T,W,H,header);


///declararea vectorului si salvarea tuturor corelatiilor celor 10 sabloane in acesta
    double *D;
    int nrd=0;
    D=(double*)malloc(10*(W-11)*(H-15)*sizeof(double));

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr0[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr1[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr2[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr3[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr4[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr5[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr6[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr7[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr8[i][j];

    for(i=0; i<(H-15); i++)
        for(j=0; j<(W-11); j++)
            D[nrd++]=corr9[i][j];

///SORTAREA detectiilor in functie de corelatie
    qsort(D,10*(W-11)*(H-15),8,sort);


///dezalocare memoriei folosite pentru partea de template matching din main

    for(i=0; i<H; i++)
        free(I[i]);
    free(I);
    for(i=0; i<H; i++)
        free(T[i]);
    free(T);
    for(i=0; i<H-15; i++)
        free(corr0[i]);
    free(corr0);
    for(i=0; i<H-15; i++)
        free(corr1[i]);
    free(corr1);
    for(i=0; i<H-15; i++)
        free(corr2[i]);
    free(corr2);
    for(i=0; i<H-15; i++)
        free(corr3[i]);
    free(corr3);
    for(i=0; i<H-15; i++)
        free(corr4[i]);
    free(corr4);
    for(i=0; i<H-15; i++)
        free(corr5[i]);
    free(corr5);
    for(i=0; i<H-15; i++)
        free(corr6[i]);
    free(corr6);
    for(i=0; i<H-15; i++)
        free(corr7[i]);
    free(corr7);
    for(i=0; i<H-15; i++)
        free(corr8[i]);
    free(corr8);
    for(i=0; i<H-15; i++)
        free(corr9[i]);
    free(corr9);
    free(D);




return 0;
}
