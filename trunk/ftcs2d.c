#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int lattice_size_x,lattice_size_y,*iup,*idown,*jup,*jdown;

#define FREE_ARG char*
#define NR_END 1

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

void setupgridneighbors()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=1;
     iup[lattice_size_x]=lattice_size_x;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

main()
{    FILE *fp1,*fp2;
     float delta,**topo,**topoold,D,duration,timestep;
     int i,j,t,nsteps;

     fp1=fopen("inputdem","r");
     fp2=fopen("outputdem","w");
     lattice_size_x=300;
     lattice_size_y=300;
     delta=10.0;     /* m */
     D=1.0;          /* m^2/kyr */
     duration=100.0; /* kyr */
     setupgridneighbors();
     topo=matrix(1,lattice_size_x,1,lattice_size_y);
     topoold=matrix(1,lattice_size_x,1,lattice_size_y);
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fscanf(fp1,"%f",&topo[i][j]);
        topoold[i][j]=topo[i][j];}
     timestep=0.5*delta*delta/(2*D);
     nsteps=(int)(duration/timestep);
     for (t=1;t<=nsteps;t++)
      {for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
          topo[i][j]+=timestep*D/(delta*delta)*(topoold[iup[i]][j]+topoold[idown[i]][j]
           +topoold[i][jup[j]]+topoold[i][jdown[j]]-4*topoold[i][j]);
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         topoold[i][j]=topo[i][j];}
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fprintf(fp2,"%f\n",topo[i][j]);
     fclose(fp1);
     fclose(fp2);
}
