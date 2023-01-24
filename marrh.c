#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>


#define  PICSIZE 256
#define  MAXMASK 100

         int    pic[PICSIZE][PICSIZE];
         double outpic1[PICSIZE][PICSIZE];
         double outpic2[PICSIZE][PICSIZE];
         int    edgeflag[PICSIZE][PICSIZE];
         double xmask[MAXMASK][MAXMASK];
         double ymask[MAXMASK][MAXMASK];
         double conv1[PICSIZE][PICSIZE];
         double conv2[PICSIZE][PICSIZE];
         double ival[256][256],maxival;


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
        double  maskval,sum1, sum2,sig,maxival,minival,maxval,ZEROTOL;
        FILE    *fo1, *fo2,*fp1, *fopen();
        char    *foobar;

        //Take the inputs and define outputs
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

        argc--; argv++;
        foobar = *argv;
        fo1=fopen(foobar,"wb");

        // argc--; argv++;
        // foobar = *argv;
        // fo2=fopen(foobar,"wb");

        //Input Sigma
        argc--; argv++;
        foobar = *argv;
        sig = atoi(foobar);
        printf("%f",sig);

        //Input Tolerance
        argc--; argv++;
        foobar = *argv;
        ZEROTOL = atoi(foobar);

        //Mask ratio aka Sigma 
        mr = (int)(sig * 3);
        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        
        initializeFiles(fo1);


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
             outpic1[i][j] = sum1;
             outpic2[i][j] = sum2;
             conv1[i][j] = sum1;
             conv2[i][j] = sum2;
          }
        }

       //Sqrt code from sobel to get magnitude, scale the output and print it out 
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             ival[i][j]=sqrt((double)((outpic1[i][j]*outpic1[i][j]) +
                                      (outpic2[i][j]*outpic2[i][j])));
             if (ival[i][j] > maxival)
                maxival = ival[i][j];

           }
        }

        //Magnitude image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             ival[i][j] = (ival[i][j] / maxival) * 255;
             fprintf(fo1,"%c",(char)((int)(ival[i][j])));
            }
          }
}
