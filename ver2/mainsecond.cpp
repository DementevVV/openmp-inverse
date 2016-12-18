#include <stdio.h>
#include <sys/time.h>
#define SIZE_NUMBER 2000

//чтение матрицы из файла
int read_file_TO_matrix(char *name,int &str,int &stb,double **&mass) {
  char buffer;
  double number;
  double ttt=10;
  int i,j;
  bool flag=false;
  mass=0;
  str=stb=i=j=0;
  unsigned int code=0;
  FILE *fo;
  fo=fopen(name,"r");
  if(!fo)
  return 0;
  while(!feof(fo)) {
    if(!fread(&buffer,sizeof(char),1,fo))
	break;
	if(buffer>='0' && buffer<='9') {
	  if(!(code>0 && code<=SIZE_NUMBER))
	    code=1;
	  else code++;
	  if(ttt==10)
	    number=number*ttt+(buffer-'0');
	  else {
	    number=number+(buffer-'0')*ttt;
		ttt/=10;
	  }
	  continue;
    }
  else if(buffer=='-') {
    flag=true;
	continue ;
  }
  else if(buffer=='.') {
    if(ttt==10) {
	  ttt=0.1;
	  continue;
	}
    else {
      for(int k=0;k<i;k++)
	    delete []mass[k]; 
	    if(mass!=0)delete []mass;
	    return 0;
    }
  }
    else if(code>0 && code<=SIZE_NUMBER) {
      if(flag)
	  number=-number;
	  if(str==0)
	  str=number;
	  else if(stb==0)
	  stb=number;
	  else if(str>0 && stb>0) {
	    if(i==0 && j==0) {
	      mass=new double *[str];
		  mass[i]=new double [stb];
	    }
	  else if(j>=stb) {
	    if(++i>=str){ j=-1;break ; }
	    mass[i]=new double [stb];
	    j=0;
	  }
	    if(j<stb)
	    mass[i][j++]=number;
	  }
    }
    ttt=10;
    flag=false;
    number=0;
    code=0;
    }
  if(str-1==i && stb==j)
  return 1;
  for(j=i;j>=0;j--)
  delete []mass[j];
  delete []mass;
  return 0;
}

//вывод матрицы в консоль
void print_matrix(int str,int stb,double **mass) {
  for(int i=0;i<str;i++) {
    for(int j=0;j<stb;j++)
	  printf("a[%d][%d]=%f \t",i,j,mass[i][j]);
	  printf("\n");
  }
  printf("\n");
}

int LU_fun_OBR(int cnt_str,double **mass,double **&LU,double **&M_obr) {
  int i,j,k;
  double sum;
  LU=new double *[cnt_str];//массив под матрицу LU
  for(i=0;i<cnt_str;i++)
    LU[i]=new double [cnt_str];

  for(i=0;i<cnt_str;i++) {
    for(j=0;j<cnt_str;j++) {
	  sum=0;
	  if(i<=j) {
	    for(k=0;k<i;k++)
		  sum+=LU[i][k]*LU[k][j];
		  //вычисл эл-ты верхней треуг матр
		  LU[i][j]=mass[i][j]-sum;
	  }
	  else {
	    for(k=0;k<j;k++)
		  sum+=LU[i][k]*LU[k][j];
		  if(LU[j][j]==0)
		    return 0;
		    //вычисл эл-ты нижней треуг матр
			LU[i][j]=(mass[i][j]-sum)/LU[j][j];
	  }
	}
}

  M_obr=new double *[cnt_str]; //созд масс для обр матр
  for(i=0;i<cnt_str;i++)
    M_obr[i]=new double [cnt_str];
  int p;
  //нахожд обр матр
  for(i=cnt_str-1;i>=0;i--) {
    for(j=cnt_str-1;j>=0;j--) {
      sum=0;
	  if(i==j) {
	    for(p=j+1;p<cnt_str;p++)
	      sum+=LU[j][p]*M_obr[p][j];
		  M_obr[j][j]=(1-sum)/LU[j][j];
	  }
	  else if(i<j) {
	    for(p=i+1;p<cnt_str;p++)
	      sum+=LU[i][p]*M_obr[p][j];
	      M_obr[i][j]=-sum/LU[i][i];
	  }
	  else {
	    for(p=j+1;p<cnt_str;p++)
	      sum+=M_obr[i][p]*LU[p][j];
		  M_obr[i][j]=-sum;
	  }
    }
  }
  return 1;
}

int main(int args,char *argv[]) {
  int str,stb;
  str=stb=0;
  double **mass;
  double **M_obr;
  double **LU;
  //чтение матрицы из файла
  if(read_file_TO_matrix(argv[1],str,stb,mass)==1 && str==stb) {
    struct timeval begin, end;
    gettimeofday(&begin, NULL);
    //разлож матр на LU и обращ матр
		if(LU_fun_OBR(str,mass,LU,M_obr)==1) {
		  gettimeofday(&end, NULL);
          //printf("Вывод LU матрицы\n");
		  //print_matrix(stb,stb,LU);//вывод LU матрицы
          //printf("Вывод обращённой матрицы\n");
		  //print_matrix(stb,stb,M_obr);//вывод обр матрицы
		  printf("%f\n", (double)((end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec) * 0.000001));
		}
		else { 
		  printf("Error solve matrix\n");
		}
    for(int i=0;i<str;i++){delete []mass[i];delete []LU[i];delete []M_obr[i];}
	delete []mass;
	delete []LU;
	delete []M_obr;
  }
  else {
    printf("Error read matrix\n");
  }
  return 0;
}