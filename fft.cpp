#include <vector>
#include <iostream>
#include <time.h>
#include <random>
#define uint unsigned int
const double dpi=3.14159265358979323846264338328;
const double dpi2=dpi*2;
std::random_device rd;
std::mt19937 mt(rd()+(uint)time(NULL));

uint bitRev(uint a)
{
  a = (a & 0x55555555) << 1 | (a & 0xAAAAAAAA) >> 1;
  a = (a & 0x33333333) << 2 | (a & 0xCCCCCCCC) >> 2;
  a = (a & 0x0F0F0F0F) << 4 | (a & 0xF0F0F0F0) >> 4;
  a = (a & 0x00FF00FF) << 8 | (a & 0xFF00FF00) >> 8;
  a = (a & 0x0000FFFF) << 16 | (a & 0xFFFF0000) >> 16;
  return a;
}

//inv=1 fft,inv=-1 ifft
void fft(std::vector<double> &a,std::vector<double> &b,int inv)
{
  uint digitN=a.size();//2のべき乗である必要がある
  double radbase;
  uint digit_level=0;
  for(uint i=1;i<digitN;i*=2)
  {
    digit_level++;
    uint blocksz=digitN/2/i;
    radbase=dpi/blocksz*inv;//事前計算できるところは全部する
    for(uint j=0;j<digitN/2;j++)
    {
      uint t2 = j%blocksz;
      uint t0 = (j/blocksz)*blocksz*2 + t2;
      uint t1 = t0+blocksz;
      double rx=a[t1];
      double ix=b[t1];
      double rm0=a[t0]-rx;
      double im0=b[t0]-ix;
      double dsin=sin(radbase*t2);
      double dcos=cos(radbase*t2);
      a[t0] += rx;
      b[t0] += ix;
      a[t1] = rm0*dcos-im0*dsin;
      b[t1] = rm0*dsin+im0*dcos;
    }
  }

  double rd=1.0;//ビット逆順して1/nを適宜かける
  if (inv==-1)rd=1.0/digitN;
  for(uint i=0;i<digitN;i++)
  {
    uint ri=bitRev(i)>>(32-digit_level);
    if (i<=ri){//<じゃなく<=じゃないとダメ
      double tmp=a[ri];
      a[ri]=a[i]*rd;
      a[i]=tmp*rd;
      tmp=b[ri];
      b[ri]=b[i]*rd;
      b[i]=tmp*rd;
    }
  }
}


//inr → outr,outi
void fftr2c(std::vector<double> &inr,std::vector<double> &outr,std::vector<double> &outi,int inv){
  int digitN=inr.size()/2;//2のべき乗である必要がある
  std::vector<double> A_(digitN),B_(digitN);
  for(int i=0;i<digitN;i++){
    A_[i]=inr[i*2];
    B_[i]=inr[i*2+1];
  }
  
  fft(A_,B_,inv);//ここが本来の1/2のサイズのfftになる

  for(int k=1;k<digitN;k++){
    double tr=A_[k]-A_[digitN-k];
    double ti=B_[k]+B_[digitN-k];
    double ur=-sin(dpi*k/digitN);
    double ui=cos(dpi*k/digitN);
    double sr=tr*ur-ti*ui;
    double si=ti*ur+tr*ui;
    outr[k]=-0.5*sr+0.5*(A_[k]+A_[digitN-k]);
    outi[k]=-0.5*si+0.5*(B_[k]-B_[digitN-k]);
  }
  outr[0]=A_[0]+B_[0];
  outr[digitN]=A_[0]-B_[0];

  for(int i=1;i<digitN;i++)
    outr[digitN+i]=outr[digitN-i];
  for(int i=1;i<digitN;i++)
    outi[digitN+i]=-outi[digitN-i];
}


//inr,ini → outr
void fftc2r(std::vector<double> inr,std::vector<double> ini,std::vector<double> &outr,int inv){
  uint digitN=inr.size()/2;//2のべき乗である必要がある
  std::vector<double> A_(digitN),B_(digitN);
  
  //inr,iniはこのように複素共役で左右対称である必要がある
  for(int i=1;i<digitN;i++)
    inr[digitN-i]=0.5*(inr[digitN-i]+inr[digitN+i]);
  for(int i=1;i<digitN;i++)
    ini[digitN-i]=0.5*(ini[digitN-i]-ini[digitN+i]);
  ini[0]=0;ini[digitN]=0;
  A_[0]=0.5*(inr[0]+inr[digitN]);
  B_[0]=0.5*(inr[0]-inr[digitN]);

  for(int k=1;k<digitN;k++){
    double tr=inr[k]-inr[digitN-k];
    double ti=ini[k]+ini[digitN-k];

    double ur=sin(dpi*k/digitN);
    double ui=cos(dpi*k/digitN);
    double sr=tr*ur-ti*ui;
    double si=ti*ur+tr*ui;

    A_[k]=0.5*sr+0.5*(inr[k]+inr[digitN-k]);
    B_[k]=0.5*si+0.5*(ini[k]-ini[digitN-k]);
  }

  fft(A_,B_,inv);//ここが本来の1/2のサイズのfftになる
  for(int i=0;i<digitN;i++){
    outr[i*2]=A_[i];
    outr[i*2+1]=B_[i];
  }
}


int main(){
  uint digitN=16;//2のべき乗である必要がある
  std::vector<double> A(digitN),B(digitN),Ac2r(digitN),Bc2r(digitN),R(digitN);
  
  for(int i=0;i<digitN;i++){
    Ac2r[i]=A[i]=(double)(mt()%65536)/256.0;//ランダム初期値生成
    Bc2r[i]=B[i]=(double)(mt()%65536)/256.0;//ランダム初期値生成
  }

  fftc2r(Ac2r,Bc2r,R,-1);
  fft(A,B,-1);
  
  std::cout<<"iFFT result"<<std::endl;
  for (auto elem : A) {//普通のiFFTの実数の結果
    std::cout<<elem<<" ";
  }
  std::cout<<std::endl;

  std::cout<<"iFFTc2r result"<<std::endl;
  for (auto elem : R) {//c2rで計算した実数の結果
    std::cout<<elem<<" ";
  }
  std::cout<<std::endl;


  double rr=0.0;
  for(int i=0;i<digitN;i++){
    rr+=abs(A[i]-R[i]);
  }
  std::cout<<"All residual"<<std::endl;
  std::cout<<rr<<std::endl;
  return 0;
}