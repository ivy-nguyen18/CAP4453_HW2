#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>


#define  PICSIZE 256
#define  MAXMASK 100

         int    pic[PICSIZE][PICSIZE];
         double outpicX[PICSIZE][PICSIZE];
         double outpicY[PICSIZE][PICSIZE];
         int    edgeflag[PICSIZE][PICSIZE];
         double xmask[MAXMASK][MAXMASK];
         double ymask[MAXMASK][MAXMASK];
         double conv_x[PICSIZE][PICSIZE];
         double conv_y[PICSIZE][PICSIZE];
         double mag[256][256],maxmag, cand[256][256], final[256][256];
         int hist[256];


void initializeFiles(FILE *f){
  fprintf(f, "P5\n");
  fprintf(f, "%d %d\n", 256, 256);  
  fprintf(f, "255\n");
}

int main(argc,argv)
int argc;
char **argv;
{
        int     i,j,p,q,s,t,mr,centx,centy;
        double  maskval,sum1, sum2,sig,maxmag,minmag,maxval, LO, HI, slope, percent, cutOff, areaOfTops;
        FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
        char    *foobar;

        //Take the inputs and define outputs
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

        argc--; argv++;
        foobar = *argv;
        fo1=fopen(foobar,"wb");

        argc--; argv++;
        foobar = *argv;
        fo2=fopen(foobar,"wb");

        argc--; argv++;
        foobar = *argv;
        fo3=fopen(foobar,"wb");

        //Input Sigma
        argc--; argv++;
        foobar = *argv;
        sig = atoi(foobar);

        //Input HI Tolerance
        // argc--; argv++;
        // foobar = *argv;
        // HI = atoi(foobar);

        //Input LO Tolerance
        // argc--; argv++;
        // foobar = *argv;
        // LO = atoi(foobar);

        //Input Percentage
        argc--; argv++;
        foobar = *argv;
        percent = atof(foobar);

        //Mask ratio aka Sigma 
        mr = (int)(sig * 3);
        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        
        initializeFiles(fo1);
        initializeFiles(fo2);
        initializeFiles(fo3);


        //Skip the heading when reading file
        for(i=0; i<3; i++){
          fscanf(fp1, "%*[^\n]\n");
        }

        //Read the picture
        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                }
        }

        //Take the first order derivative - xmask and ymask
        for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
              
              (xmask[p+centy][q+centx]) = (q*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
              (ymask[p+centy][q+centx]) = (p*(exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
           }
        }

        //Convolution
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum1 = 0;
             sum2 = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum1 += pic[i+p][j+q] * xmask[p+centy][q+centx];
                   sum2 += pic[i+p][j+q] * ymask[p+centy][q+centx];

                }
             }
             outpicX[i][j] = sum1;
             outpicY[i][j] = sum2;
             conv_x[i][j] = sum1;
             conv_y[i][j] = sum2;
          }
        }

       //Sqrt code from sobel to get magnitude, scale the output and print it out 
        maxmag = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             mag[i][j]=sqrt((double)((outpicX[i][j]*outpicX[i][j]) +
                                      (outpicY[i][j]*outpicY[i][j])));
             if (mag[i][j] > maxmag)
                maxmag = mag[i][j];

           }
        }

        //Magnitude image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             mag[i][j] = (mag[i][j] / maxmag) * 255;
             fprintf(fo1,"%c",(char)((int)(mag[i][j])));
            }
          }

        //Indicate Peaks
        for(i=mr;i<256-mr;i++){
            for(j=mr;j<256-mr;j++){

                if((conv_x[i][j]) == 0.0) {
                    conv_x[i][j] = .00001;
                }
                
                slope = conv_y[i][j]/conv_x[i][j];
                if( (slope <= .4142)&&(slope > -.4142)){
                    
                    if((mag[i][j] > mag[i][j-1])&&(mag[i][j] > mag[i][j+1])){
                        cand[i][j] = 255;
                    }
                }
                else if( (slope <= 2.4142)&&(slope > .4142)){
                 if((mag[i][j] > mag[i-1][j-1])&&(mag[i][j] > mag[i+1][j+1])){
                     cand[i][j] = 255;
                   }
                }
                else if( (slope <= -.4142)&&(slope > -2.4142)){
                 if((mag[i][j] > mag[i+1][j-1])&&(mag[i][j] > mag[i-1][j+1])){
                     cand[i][j] = 255;
                   }
                }else{
                 if((mag[i][j] > mag[i-1][j])&&(mag[i][j] > mag[i+1][j])){
                     cand[i][j] = 255;
                  }
                }
            }
        }

        //Peaks image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             fprintf(fo2,"%c",(char)((int)(cand[i][j])));
            }
        }

        //Use Percentage as Hi and Lo
        for (i=mr;i<256-mr;i++)
          { for (j=mr;j<256-mr;j++)
            {
                hist[(int)mag[i][j]]+=1;
            }
          }

        cutOff = percent*256*256;
        for(i=256; i>1; i--){
            areaOfTops += hist[i];
            if(areaOfTops > cutOff) break;
        }

        HI = i;
        LO = .35*HI;
        printf("HI THRESHOLD %f", HI);

        //Find final peaks
        for(i=mr;i<256-mr;i++){
            for(j=mr;j<256-mr;j++){

                if(cand[i][j] == 255){
                    if(mag[i][j] > HI){
                        cand[i][j] = 0;
                        final[i][j] = 255;
                    }
                    else if(mag[i][j] < LO){
                        cand[i][j] = 0;
                        final[i][j] = 0;
                    }
                }
            }
        }

        int moretodo = 1;
        while(moretodo == 1){
            moretodo = 0;
            for (i=mr;i<256-mr;i++)
            { for (j=mr;j<256-mr;j++)
                {
                    if (cand[i][j] == 255){
                        for(p=-1; p<=1; p++){
                            for(q=-1; q<=1; q++){
                                if(final[i+p][j+q] == 255){
                                    cand[i][j] = 0;
                                    final[i][j] = 255;
                                    moretodo = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        //Final image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             fprintf(fo3,"%c",(char)((int)(final[i][j])));
            }
          }
}
