/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    unsigned int i,j,z;
    ip_mat *a=(ip_mat*)malloc(sizeof(ip_mat)); /*alloco lo spazio necessario per una ip_mat*/
    a->w=w;
    a->h=h;
    a->k=k;
    a->stat=(stats*)malloc(sizeof(stats)*k);
    a->data=(float***)malloc(sizeof(float*)*h);
    for(i=0;i<h;i++){
        a->data[i]=(float**)malloc(sizeof(float*)*w);
        for(j=0;j<w;j++){
            a->data[i][j]=(float*)malloc(sizeof(float*)*k); /*qui*/
            for(z=0;z<k;z++){
                a->data[i][j][z]=v;
            }
        }
    }
    return a;
}
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    ip_mat * out=NULL;
    unsigned int i,j,z;
    /*scorro e unisco le due immagini per la dimensione prefissata*/
    if(dimensione==0){
        out=ip_mat_create((a->h)+(b->h),a->w,a->k,(float)0.0);
        for(i=0;i<(a->h+b->h);i++){
            for(j=0;j<(a->w);j++){
                for(z=0;z<(a->k);z++){
                    if(i<(a->h)) {
                        out->data[i][j][z] = a->data[i][j][z];
                    }else{
                        out->data[i][j][z] = b->data[i-(a->h)][j][z];
                    }
                }
            }
        }
    }else if(dimensione==1){
        out=ip_mat_create(a->h,a->w+b->w,a->k,(float)0.0);
        for(i=0;i<(a->h);i++){
            for(j=0;j<(a->w+b->w);j++){
                for(z=0;z<(a->k);z++){
                    if(j<(a->w)) {
                        out->data[i][j][z] = a->data[i][j][z];
                    }else{
                        out->data[i][j][z] = b->data[i][j-(a->w)][z];
                    }
                }
            }
        }
    }else if(dimensione==2){
        out=ip_mat_create(a->h,a->w,a->k+b->k ,(float)0.0);
        for(i=0;i<(a->h);i++){
            for(j=0;j<(a->w);j++){
                for(z=0;z<(a->k+b->k);z++){
                    if(z<(a->k)) {
                        out->data[i][j][z] = a->data[i][j][z];
                    }else{
                        out->data[i][j][z] = b->data[i][j][z-(a->k)];
                    }
                }
            }
        }
    }
    return out;
}
ip_mat * ip_mat_copy(ip_mat * in){
    ip_mat *strut=ip_mat_create(in->h, in->w, in->k, (float)0.0);
        strut->stat=in->stat;
        strut->data=in->data;
    return strut;

}
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    /*faccio la somma dei valori di due ipmat che ho in input*/
    ip_mat *s=ip_mat_create(a->h, a->w, a->k, (float)0.0);
    unsigned int i,j,l;
    for(i=0;i<(a->h) && i<(b->h);i++){
        for(j=0;j<(a->w) && j<(b->w);j++){
            for(l=0;l<(a->k) && l<(b->k);l++){
                s->data[i][j][l]=((a->data[i][j][l])+(b->data[i][j][l]));
            }
        }
    }
    compute_stats(s);
    return s;
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    ip_mat *sub =ip_mat_create(a->h, a->w, a->k, (float)0.0);
    unsigned int i,j,l;
    for(i=0;i<(a->h) && i<(b->h);i++){
        for(j=0;j<(a->w) && j<(b->w);j++){
            for(l=0;l<(a->k) && l<(b->k);l++){
                sub->data[i][j][l]=((a->data[i][j][l])-(b->data[i][j][l]));
            }
        }
    }

    compute_stats(sub);
    return sub;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    ip_mat *media=ip_mat_create(a->h, a->w, a->k, (float)0.0);
    unsigned int i,j,l;
    int acc=0;
    for(i=0;i<(a->h) && i<(b->h);i++){
        for(j=0;j<(a->w) && j<(b->w);j++){
            for(l=0;l<(a->k) && l<(b->k);l++){
                media->data[i][j][l]=((a->data[i][j][l])+(b->data[i][j][l]));
                acc++;
            }
        }
    }
    for(i=0;i<(a->h) && i<(b->h);i++){
        for(j=0;j<(a->w) && j<(b->w);j++){
            for(l=0;l<(a->k) && l<(b->k);l++){
                media->data[i][j][l]=media->data[i][j][l]/(float)acc;
            }
        }
    }
    compute_stats(media);
    return media;
}
void ip_mat_free(ip_mat *a){
    unsigned int i,j,z,x;
    if(a!=NULL) {   /*se si può libero la memoria*/
        free(a->stat);  /*libero array di stat*/
        for (i = 0; i < (a->h); i++) {
            for (j = 0; j < (a->w); j++) {
                free(a->data[i][j]);    /*libero tutti gli z*/
            }
            free(a->data[i]);   /*libero tutti gli i*/
        }
        free(a->data);  /*libero la matrice*/
        free(a);    /*libero tutto */
    }
}

void compute_stats(ip_mat * t){
    unsigned int i, j, z, q, cont=0, a;
    float *min=(float*)malloc((t->k) * sizeof(float));
    float *max=(float*)malloc((t->k) * sizeof(float));
    float *mean=(float*)malloc((t->k) * sizeof(float));
    for(q=0;q<(t->k);q++){
        min[q]=t->data[0][0][q];
        max[q]=t->data[0][0][q];
    }
    for(i=0;i<(t->h);i++){
        for(j=0;j<(t->w);j++){
            for(z=0;z<(t->k);z++){
                if(min[z]>(t->data[i][j][z])){
                    min[z]=t->data[i][j][z];
                }
                if(max[z]<(t->data[i][j][z])){
                    max[z]=t->data[i][j][z];
                }
                if(z==0){
                    cont++;
                }
                mean[z]=mean[z]+t->data[i][j][z];
            }
        }
    }
    for(q=0;q<(t->k);q++){
        mean[q]=mean[q]/(float)cont;
    }
    for(a=0;a<(t->k);a++){
        t->stat[a].min=min[a];
        t->stat[a].max=max[a];
        t->stat[a].mean=mean[a];
    }
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    unsigned int i,j, row_dim=row_end-row_start+1, col_dim=col_end-col_start+1;
    unsigned int salvaColStart=col_start;   /*mi servirà nel ciclo per rimettere la variabile alla dimensione originale*/
    ip_mat *sub=ip_mat_create(row_dim, col_dim, t->k, (float)0.0);
    for(i=0;i<row_dim;i++, ++row_start){
        col_start=salvaColStart;
        for(j=0;j<col_dim;j++, ++col_start){
            sub->data[i][j]=t->data[row_start][col_start];
        }
    }
    compute_stats(sub);
    return sub;
}

void ip_mat_init_random(ip_mat * t, float mean, float var){
    unsigned int i,j,k;
    float sigma=(float)sqrt((double)var);/*calcolo il sigma con la radice quadrata della varianza*/
    t=ip_mat_create(t->h, t->w, t->k, (float)0.0);
    for(i=0;i<(t->h);i++){
        for(j=0;j<(t->w);j++){
            for(k=0;k<(t->k);k++){
                t->data[i][j][k]=get_normal_random(mean, sigma);    /*chiamo la funzione e riempio la mtrice di numeri casuali*/
            }
        }
    }
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    ip_mat *tensore =ip_mat_create(a->h, a->w, a->k, (float)0.0);   /*variabile di supporto*/
    unsigned int i,j,k;
    for(i=0;i<(a->h);i++){
        for(j=0;j<(a->w);j++){
            for(k=0;k<(a->k);k++){
                tensore->data[i][j][k]=(a->data[i][j][k])+c;
            }
        }
    }
    return tensore;
}
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    ip_mat *tensore =ip_mat_create(a->h, a->w, a->k, (float)0.0);
    unsigned int i,j,k;
    for(i=0;i<(a->h);i++){
        for(j=0;j<(a->w);j++){
            for(k=0;k<(a->k);k++){
                tensore->data[i][j][k]=(a->data[i][j][k])*c;
            }
        }
    }
    return tensore;
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    if(in!=NULL) {
        ip_mat *s = ip_mat_create(in->h, in->w, in->k, (float) 0.0);
        unsigned int i, j, k;
        float sum = 0, media = 0;
        for (i = 0; i < (in->h); i++) {
            for (j = 0; j < (in->w); j++) {
                for (k = 0; k < (in->k); k++) {
                    sum = sum + (in->data[i][j][k]);
                    /*faccio una variabile che contiene il risultato della somma di tutti i canali
                    *presenti in quel pixel*/
                }
                media = sum / (float) k;    /*faccio la media quindi somma diviso il numero di canali*/
                for (k = 0; k < (in->k); k++) {
                    s->data[i][j][k] = media;   /*metto la media dentro la ip_mat*/
                }
                sum = 0;
            }
        }
        return s;
    }else{
        exit(1);
    }


}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    if(a->h == b->h && a->w == b->w && a->k == b->k) {
        /*le due immagini devono avere la stessa dimensione assolutamente*/
        ip_mat *new = ip_mat_create(a->h, a->w, a->k, (float) 0.0);
        unsigned int i, j, k;
        for (i = 0; i < (a->h); ++i) {
            for (j = 0; j < (a->w); ++j) {
                for (k = 0; k < (a->k); ++k) {
                    new->data[i][j][k] = ((a->data[i][j][k]) * alpha) + ((b->data[i][j][k]) * (1 - alpha));
                    /*tramite la formula metto nella mia nuova immagine la somma della prima*alpha + la seconda
                     * moltiplicata *(1-alpha)*/
                }
            }
        }
        return new;
    }else{
        exit(1);
    }
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    if(a!=NULL) {
        ip_mat *new = ip_mat_create(a->h, a->w, a->k, (float) 0.0);
        unsigned int i, j, k;
        for (i = 0; i < (a->h); i++) {
            for (j = 0; j < (a->w); j++) {
                for (k = 0; k < (a->k); k++) {
                    new->data[i][j][k] = (a->data[i][j][k]) + bright;
                    /*semplicemente aggiungo ad ogni pixel una determinata luminosità*/
                }
            }
        }
        return new;
    }else{
        exit(1);
    }
}
ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    if(a!=NULL) {
        ip_mat *out = ip_mat_copy(a);   /*me la copio in un altra variabile così non modifico la originale*/
        float mean = (float) 0.0, var = (float) 0.1;    /*setto media e varianza*/
        ip_mat_init_random(out, mean, var); /*inizializzo la ipmat con numeri random*/
        ip_mat_mul_scalar(out, amount); /*moltiplico la matrice per il numero di rumore che voglio*/
        out = ip_mat_sum(a, out);   /*sommo la immagine iniziale con quella che ho generato casualmente e metto in out il risultato*/
        return out;
    }else{
        exit(1);
    }
}
float moltiplica_per_filtro(ip_mat *sott, ip_mat *filtro, int k){
    /*devo moltiplicare due matrici fra loro una è la sottomatrice che ha dimensioni del filtro e l'altra è il filtro*/
    unsigned int i,j;
    float risultato=(float)0.0;
    for(i=0;i<(sott->h);i++){
        for(j=0;j<(sott->w);j++){
                risultato+=((sott->data[i][j][k])*(filtro->data[i][j][k]));
        }
    }
    /*quello che ottengo è un numero solo che è la somma di tutte le moltiplicazioni che ho fatto fra la matrice 1 e il filtro*/
    return risultato;
}
ip_mat * crea_sotto_matrice_da_indice_di_dimensione(unsigned int h, unsigned int w, ip_mat *original, unsigned int dimFiltroh, unsigned int dimFiltrow,unsigned int dimFiltrok){
    ip_mat *temp=ip_mat_subset(original,h,h+dimFiltroh-1,w,w+dimFiltrow-1);
    /*chiamo subset e faccio una sottomatrice della dimensione del filtro così la posso moltiplicare con il filtro*/
    return temp;
}
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    ip_mat *out = ip_mat_create(a->h, a->w, a->k, (float)0.0);
    ip_mat *a_con_padding = ip_mat_padding(a, (((f->h)-1)/2), (((f->w)-1)/2));
    /*aggiungo all'immaginina iniziale un bordo bianco e ci metto l'immagine dentro, il bordo dipende dalla dimensione del filtro
     *non dall'immagine originale*/
    unsigned int i,j,k;
    for(i=0;i<(a_con_padding->h)-((f->h)-1);i++) {
        for (j = 0; j < (a_con_padding->w) - ((f->w) - 1); j++) {
            for (k = 0; k < (a_con_padding->k); k++) {
                ip_mat *sottomatrice = crea_sotto_matrice_da_indice_di_dimensione(i, j, a_con_padding, f->h, f->w, f->k);
                /*creo una sottomatrice a seconda della posizione in cui sono*/
                out->data[i][j][k] = moltiplica_per_filtro(sottomatrice, f, k);
                /*moltiplico la sottomatrice per il filtro e la metto nella immagine risultante*/
                /*quindi l'immagine risultante avrà la stessa dimensione di quella originale*/
            }
        }
    }
    /*scompongo la matrice in tante sottomatrici della dimensione di f e moltiplico per i valori di f
     * e salvo il singolo valore nella matrice out*/
    ip_mat_free(a_con_padding);/*libero la ip_mat col bordino che avevo usato solo per moltiplicarla con il filtro*/
    return out;
}

ip_mat * create_edge_filter(){
    int i, j, z;
    ip_mat *new=ip_mat_create(3,3,3 , (float)0.0);
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            for(z=0; z<3; z++) {
                if (i == 1 && j ==1 ) {
                    new->data[i][j][z] = (float) 8.0;
                } else {
                    new->data[i][j][z] = (float) -1.0;
                }

            }
        }
    }
    return new;
}

ip_mat * create_emboss_filter(){
    int i, j, z;
    ip_mat *new=ip_mat_create(3,3,3 , (float)0.0);

    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            for(z=0; z<3; z++) {
                if ((i==0 && j == 1)||(i==1 && j==0)) {
                    new->data[i][j][z] = (float) -1.0;
                } else if ((i==0 && j == 2)||(i==2 && j==0)) {
                    new->data[i][j][z] = (float) 0.0;
                } else if (i==0 && j == 0) {
                    new->data[i][j][z] = (float) -2.0;
                } else if (i==2 && j == 2) {
                    new->data[i][j][z] = (float) 2.0;
                } else {
                    new->data[i][j][z] = (float) 1.0;
                }
            }
        }
    }
    return new;
}
/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
    unsigned int i,j,k;
    for(i=0;i<(t->h);i++){
        for(j=0;j<(t->w);j++){
            for(k=0;k<(t->k);k++){
                if(t->data[i][j][k]<low){
                    t->data[i][j][k]=low;
                }else if(t->data[i][j][k]>high){
                    t->data[i][j][k]=high;
                }
            }
        }
    }
}

ip_mat * create_sharpen_filter(){
    int i,j,z;
    ip_mat *new=ip_mat_create(3,3,3 , (float)0.0);
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            for(z=0; z<3; z++) {
                if ((( i == 0 || i == 2 ) && j == 0 )||(( i == 0 || i== 2 ) && j == 2 )) {
                    new->data[i][j][z] = (float) 0.0;
                } else {
                    if((( i == 0 || i == 2 ) && j == 1 ) || ( i == 1 && ( j == 0 || j == 2 ))) {
                        new->data[i][j][z] = (float) -1.0;
                    }
                    else {
                        new->data[i][j][z] = (float) 5.0;
                    }
                }

            }
        }
    }
    return new;
}

void rescale(ip_mat * t, float new_max) {
    unsigned int i,j,k;
    compute_stats(t);
    for(i=0;i<(t->h);i++){
        for(j=0;j<(t->w);j++){
            for(k=0;k<(t->k);k++){
                t->data[i][j][k]=(((t->data[i][j][k])-(t->stat[k].min))/((t->stat[k].max)-(t->stat[k].min)))*new_max;
            }
        }
    }
}
/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
 * */
ip_mat * ip_mat_padding(ip_mat * a, int pad_h, int pad_w){
    unsigned int toth=(a->h)+(2*pad_h);
    unsigned int totw=(a->w)+(2*pad_w);
    /*calcolo la dimensione totale dell'immagine compreso il bordino che si trova sia a destra
     * che a sinistra sia in alto che in basso quindi devo moltiplicare *2 sia pad_h sia pad_w*/
    ip_mat *new=ip_mat_create(toth,totw,a->k,(float)0.0);/*creo ip_mat della nuova dimensione compresa di bordo*/
    unsigned int i,j,k,k1;
    unsigned int padh= pad_h, padw=pad_w;
    /*inizializzato toth e totw altrimenti nella condizione di if mi dava errore*/
    for(i=0;i<(toth-padh);i++){
        for(j=0;j<(totw-padw);j++){
            if(i<padh || i>(toth-padh) || j<padw || j>(totw-padw)){
                /*se sono in un posto del bordino allora ci metto 0*/
                for(k=0;k<(new->k);k++){
                    new->data[i][j][k]=(float)0.0;
                }
            }else{
                /*altrimenti ci devo mettere l'immagine che trovo con l'indice meno il bordino*/
                for(k1=0;k1<(a->k);k1++) {
                    float val =a->data[i-padh][j-padw][k1];
                    new->data[i][j][k1] = val;
                }
            }
        }
    }
    return new;
}

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma){
    ip_mat *new=ip_mat_create(h,w,k,(float)0.0);    /*creo ip_mat nuova*/
    unsigned int i,j,z;
    unsigned int centerh=h/2;   /*trovo il centro*/
    unsigned int centerw=w/2;
    unsigned int a=0;
    float *somma=(float*)(malloc(k*sizeof(float))); /*creo un'array che mi serve per avere la somma per ogni canale*/
    for(a=0;a<k;a++){
        somma[a]=(float)0.0;    /*inizializzo l'array a 0*/
    }
    for(i=0;i<(new->h);i++) {
        for (j = 0; j < (new->w); j++) {
            for (z = 0; z < (new->k); z++) {
                unsigned int distanzaCentroi=centerh-i; /*trovo la distanza dal centro*/
                unsigned int distanzaCentroj=centerw-j;
                float sigmaquad=sigma*sigma;    /*faccio sigma al quadrato*/
                float primaParte = (float)(1/(2*PI*sigmaquad));
                float secondaParte = (float)((((distanzaCentroi*distanzaCentroi)+(distanzaCentroj*distanzaCentroj))/(2*sigmaquad))*-1);
                float numero=(float)(primaParte*((exp(secondaParte)))); /*applico la formula che ho diviso in due parti per comodità*/
                somma[z]+=numero;   /*faccio la somma di tutti i valori per ogni canale mi servirà dopo per la media*/
                new->data[i][j][z]=numero; /*il numero che calcolo lo metto dentro la matrice e alla fine farò la media*/
            }
        }
    }
    for(i=0;i<(new->h);i++) {
        for (j = 0; j < (new->w); j++) {
            for (z = 0; z < (new->k); z++) {
                new->data[i][j][z]=new->data[i][j][z]/somma[z]; /*metto la media dentro ogni pixel*/
            }
        }
    }
    return new;
}
/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(unsigned int h, unsigned int w, unsigned int k){
    float c = (float)1.0/(float)(w*h);  /*riempio tutta la matrice di 1/(w*h)*/
    ip_mat *out=ip_mat_create(h,w,k,c);
    return out;
}
