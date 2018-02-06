#ifndef _tensor2_cm_
#define _tensor2_cm_
#include <algorithm>
#include <array>
#include <iostream>
#include <math.h>

template <class Scal>
class coordinate  {
public :
  std::array<Scal, 3> data;
  coordinate & operator  +=( const coordinate &rhs ){
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    return *this;
  }
  coordinate & operator  -=( const coordinate &rhs ){
    data[0] -= rhs.data[0];
    data[1] -= rhs.data[1];
    data[2] -= rhs.data[2];
    return *this;
  }
  Scal & operator[](size_t i){return data[i];}
  const Scal & operator[](size_t i) const {return data[i];}
  
};

template <class Scal>
coordinate<Scal > operator+( const coordinate<Scal> & A, const coordinate<Scal> & B ){
  return coordinate<Scal >( {std::array<Scal ,3>{A[0]+B[0],A[1]+B[1],A[2]+B[2]}} ) ;
}

template <class Scal>
coordinate<Scal> operator-( const coordinate<Scal> & A, const coordinate<Scal> & B ){
  return{A[0]-B[0],A[1]-B[1],A[2]-B[2]} ;
}

template <class Scal>
class tensor2_9cm{
  // template < class T> tensor2_9rm( const T &in):data(std::forward(in)){}
public:
  typedef Scal value_type;
  tensor2_9cm() =default;
  // tensor2_9cm( std::array< Scal, 9> &&_data):data(std::move(_data)){}
  //Scal &operator()(size_t i, size_t j ) {return data[i+3*j];}
  //const Scal &operator()(size_t i, size_t j ) const {return data[i+3*j];}
  tensor2_9cm & operator*=(Scal alpha){
    //std::for_each( data.begin(), data.end(), [&alpha](Scal &a){a*=alpha;}); 
    for(auto & a :data) {a*=alpha;}
    /*data[0]*=alpha;
    data[1]*=alpha;
    data[2]*=alpha;
    data[3]*=alpha;
    data[4]*=alpha;
    data[5]*=alpha;
    data[6]*=alpha;
    data[7]*=alpha;
    data[8]*=alpha;
    */
    return *this;
   }
  tensor2_9cm & operator/=(Scal alpha){
    for(auto & a :data) {a/=alpha;}
    return *this;
  }
  tensor2_9cm<Scal> & operator+=(const  tensor2_9cm<Scal > &A){
    std::transform(A.data.begin(), A.data.end(), data.begin(), data.begin(),[](const Scal &a, const Scal &b){ return a+b; }  );
    return *this;
  }
  tensor2_9cm<Scal> & operator-=(const  tensor2_9cm<Scal > &A){
    std::transform(A.data.begin(), A.data.end(), data.begin(), data.begin(),[](const Scal &a, const Scal &b){ return b-a; }  );
    return *this;
  }
public :
  std::array< Scal , 9 > data;
};

template <class Scal>
Scal det(const tensor2_9cm<Scal> & _A ){
  /*return in.data[0][0] * (in.data[1][1] * in.data[2][2] - in.data[1][2] * in.data[2][1]) -
    in.data[0][1] * (in.data[1][0] * in.data[2][2] - in.data[1][2] * in.data[2][0]) +
    in.data[0][2] * (in.data[1][0] * in.data[2][1] - in.data[1][1] * in.data[2][0]);
  */
  const auto &A = _A.data;
  return    A[0] * (A[4] * A[8] - A[5] * A[7]) 
          - A[3] * (A[1] * A[8] - A[2] * A[7]) 
          + A[6] * (A[1] * A[5] - A[2] * A[4]);
}


template <class Scal>
tensor2_9cm<Scal> invert(const tensor2_9cm<Scal> & _A){
  const Scal detA = det(_A);
  const auto &A = _A.data;
  /*if (detin == 0.0){
    throw -1;
    }*/
  return tensor2_9cm<Scal>
    { ( A[4] * A[8] - A[5] * A[7] )/detA,
	( A[2] * A[7] - A[1] * A[8] )/detA,
	( A[1] * A[5] - A[2] * A[4] )/detA,
	( A[6] * A[5] - A[3] * A[8] )/detA,
	( A[0] * A[8] - A[2] * A[6] )/detA,
	( A[2] * A[3] - A[0] * A[5] )/detA,
	( A[3] * A[7] - A[4] * A[6] )/detA,
	( A[1] * A[6] - A[0] * A[7] )/detA,
	( A[0] * A[4] - A[1] * A[3] )/detA};
}
template <class Scal>
tensor2_9cm<Scal> invert_transpose(const tensor2_9cm<Scal> & _A){
  const Scal detA = det(_A);
  const auto &A = _A.data;
  /*if (detA == 0.0){
    throw -1;
    }*/
  return tensor2_9cm<Scal>
    {
      /*  inv00*/   ( A[4] * A[8] - A[5] * A[7] )/detA,
	/*inv01*/	( A[6] * A[5] - A[3] * A[8] )/detA,
	/*inv02*/	( A[3] * A[7] - A[4] * A[6] )/detA,
	/*inv10*/ ( A[2] * A[7] - A[1] * A[8] )/detA,
	/*inv11*/	( A[0] * A[8] - A[2] * A[6] )/detA,
	/*inv12*/	( A[1] * A[6] - A[0] * A[7] )/detA,
	/*inv20*/	( A[1] * A[5] - A[2] * A[4] )/detA,
	/*inv21*/	( A[2] * A[3] - A[0] * A[5] )/detA,
	/*inv22*/	( A[0] * A[4] - A[1] * A[3] )/detA};
}

template <class Scal>
tensor2_9cm<Scal> operator*( const tensor2_9cm<Scal> & _A,  const tensor2_9cm<Scal> & _B){
  const auto & A = _A.data;
  const auto & B = _B.data;
  
  return tensor2_9cm<Scal>{
    A[0] * B[0] + A[3] *B[1] + A[6]*B[2],
      A[1] * B[0] + A[4] *B[1] + A[7]*B[2],
      A[2] * B[0] + A[5] *B[1] + A[8]*B[2],
      A[0] * B[3] + A[3] *B[4] + A[6]*B[5],
      A[1] * B[3] + A[4] *B[4] + A[7]*B[5],
      A[2] * B[3] + A[5] *B[4] + A[8]*B[5],
      A[0] * B[6] + A[3] *B[7] + A[6]*B[8],
      A[1] * B[6] + A[4] *B[7] + A[7]*B[8],
      A[2] * B[6] + A[5] *B[7] + A[8]*B[8] };
}

template <class Scal>
coordinate<Scal> operator*( const tensor2_9cm<Scal> & _A,  const coordinate<Scal> & _b){
  const auto & A = _A.data;
  const auto & b = _b.data;
  
  return coordinate<Scal>{
    A[0] * b[0] + A[3] *b[1] + A[6]*b[2],
      A[1] * b[0] + A[4] *b[1] + A[7]*b[2],
      A[2] * b[0] + A[5] *b[1] + A[8]*b[2]
      };
}

template <class Scal >
tensor2_9cm<Scal > Green_Lagrange ( const   tensor2_9cm<Scal>  &F){
  const Scal E00 = 0.5*(F.data[0]*F.data[0] + F.data[1]*F.data[1] +F.data[2]*F.data[2] -1.);
  const Scal E10 = 0.5*(F.data[3]*F.data[0] + F.data[4]*F.data[1] + F.data[5]*F.data[2]);
  const Scal E20 = 0.5*(F.data[6]*F.data[0] + F.data[7]*F.data[1] + F.data[8]*F.data[2]);
  const Scal E11 = 0.5*(F.data[3]*F.data[3] + F.data[4]*F.data[4] + F.data[5]*F.data[5] -1.);
  const Scal E21 = 0.5*(F.data[6]*F.data[3] + F.data[7]*F.data[4] + F.data[8]*F.data[5]);
  const Scal E22 = 0.5*(F.data[6]*F.data[6] + F.data[7]*F.data[7] + F.data[8]*F.data[8]-1.);
  /*
  const Scal E00 = 0.5*(F(0,0)*F(0,0) + F(1,0)*F(1,0)+ F(2,0)*F(2,0) -1.);
  const Scal E10 = 0.5*(F(0,1)*F(0,0) + F(1,1)*F(1,0)+ F(2,1)*F(2,0) );
  const Scal E20 = 0.5*(F(0,2)*F(0,0) + F(1,2)*F(1,0)+ F(2,2)*F(2,0) );
  //const auto E01 = 0.5*(F(0,0)*F(0,1) + F(1,0)*F(1,1)+ F(2,0)*F(2,1) -1.);
  const Scal E11 = 0.5*(F(0,1)*F(0,1) + F(1,1)*F(1,1)+ F(2,1)*F(2,1) -1.);
  const Scal E21 = 0.5*(F(0,2)*F(0,1) + F(1,2)*F(1,1)+ F(2,2)*F(2,1) );
  // const auto E02 = 0.5*(F(0,0)*F(0,2) + F(1,0)*F(1,2)+ F(2,0)*F(2,2) -1.);
  //const auto E12 = 0.5*(F(0,1)*F(0,2) + F(1,1)*F(1,2)+ F(2,1)*F(2,2) -1.);
  const Scal E22 = 0.5*(F(0,2)*F(0,2) + F(1,2)*F(1,2)+ F(2,2)*F(2,2) -1.);
  //return tensor2_9cm< double> (std::array<Scal, 9> {E00, E10, E20, E01, E11, E21, E02, E12, E22 });
  //  return tensor2_9cm< Scal> (std::array<Scal, 9> {E00, E10, E20, E10, E11, E21, E20, E21, E22 });
  */
  return tensor2_9cm< Scal> {E00, E10, E20, E10, E11, E21, E20, E21, E22 };
  
  /*
  const auto E00 = 0.5*(F(0,0)*F(0,0) + F(1,0)*F(1,0)+ F(2,0)*F(2,0) -1.);
  const auto E10 = 0.5*(F(0,1)*F(0,0) + F(1,1)*F(1,0)+ F(2,1)*F(2,0) -1.);
  const auto E20 = 0.5*(F(0,2)*F(0,0) + F(1,2)*F(1,0)+ F(2,2)*F(2,0) -1.);
  const auto E01 = 0.5*(F(0,0)*F(0,1) + F(1,0)*F(1,1)+ F(2,0)*F(2,1) -1.);
  const auto E11 = 0.5*(F(0,1)*F(0,1) + F(1,1)*F(1,1)+ F(2,1)*F(2,1) -1.);
  const auto E21 = 0.5*(F(0,2)*F(0,1) + F(1,2)*F(1,1)+ F(2,2)*F(2,1) -1.);
  const auto E02 = 0.5*(F(0,0)*F(0,2) + F(1,0)*F(1,2)+ F(2,0)*F(2,2) -1.);
  const auto E12 = 0.5*(F(0,1)*F(0,2) + F(1,1)*F(1,2)+ F(2,1)*F(2,2) -1.);
  const auto E22 = 0.5*(F(0,2)*F(0,2) + F(1,2)*F(1,2)+ F(2,2)*F(2,2) -1.);
  return tensor2_9cm< Scal> (std::array<Scal, 9> {E00, E10, E20, E01, E11, E21, E02, E12, E22 });
  */
}

template <class Scal >
tensor2_9cm<Scal > Right_Cauchy_Green ( const   tensor2_9cm<Scal>  &F){
  const Scal C00 = F.data[0]*F.data[0] + F.data[1]*F.data[1] +F.data[2]*F.data[2];
  const Scal C10 = F.data[3]*F.data[0] + F.data[4]*F.data[1] + F.data[5]*F.data[2];
  const Scal C20 = F.data[6]*F.data[0] + F.data[7]*F.data[1] + F.data[8]*F.data[2];
  const Scal C11 = F.data[3]*F.data[3] + F.data[4]*F.data[4] + F.data[5]*F.data[5];
  const Scal C21 = F.data[6]*F.data[3] + F.data[7]*F.data[4] + F.data[8]*F.data[5];
  const Scal C22 = F.data[6]*F.data[6] + F.data[7]*F.data[7] + F.data[8]*F.data[8];
  return tensor2_9cm< Scal> {C00, C10, C20, C10, C11, C21, C20, C21, C22 };
}



template <class Scal>
Scal trace( const   tensor2_9cm<Scal>  &A ){
  return A.data[0] + A.data[4] + A.data[8];
}
 

template <class Scal>
// second invariant = trace ATA
Scal I2( const   tensor2_9cm<Scal>  &_A ){
  //ATA00 = A00*A00+ A10*A10 + A20*A20
  //ATA11 = A01*A01+ A11*A11 + A21*A21
  //ATA22 = A02*A02+ A12*A12 + A22*A22
  const auto & A = _A.data;
  return A[0]*A[0] + A[1]*A[1]+ A[2]*A[2]
    + A[3]*A[3] + A[4]*A[4] + A[5]*A[5]
    + A[6]*A[6] + A[7]*A[7] + A[8]*A[8];
}

template <class Scal >
Scal  norm2( const   tensor2_9cm<Scal>  &A){
  return sqrt(I2(A));
}

template <class Scal >
Scal contract2( const   tensor2_9cm<Scal>  &_A, const   tensor2_9cm<Scal>  &_B){
  const auto & A = _A.data;
  const auto & B = _B.data;
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + A[4]*B[4] + A[5]*B[5] + A[6]*B[6] + A[7]*B[7] + A[8]*B[8];  
}



/*
template <class Scal >
tensor2_9cm<Scal > operator* ( const   tensor2_9cm<Scal>  &A, const   tensor2_9cm<Scal>  &B){
  const auto C00 = A.data[0]*B.data[0] + A.data[3]*B.data[1] + A.data[6]*B.data[2];
  const auto C10 = A.data[1]*B.data[0] + A.data[4]*B.data[1] + A.data[7]*B.data[2];
  const auto C20 = A.data[2]*B.data[0] + A.data[5]*B.data[1] + A.data[8]*B.data[2];
  const auto C01 = A.data[0]*B.data[3] + A.data[3]*B.data[4] + A.data[6]*B.data[5];
  const auto C11 = A.data[1]*B.data[3] + A.data[4]*B.data[4] + A.data[7]*B.data[5];
  const auto C21 = A.data[2]*B.data[3] + A.data[5]*B.data[4] + A.data[8]*B.data[5];
  const auto C02 = A.data[0]*B.data[6] + A.data[3]*B.data[7] + A.data[6]*B.data[8];
  const auto C12 = A.data[1]*B.data[6] + A.data[4]*B.data[7] + A.data[7]*B.data[8];
  const auto C22 = A.data[2]*B.data[6] + A.data[5]*B.data[7] + A.data[8]*B.data[8];
  
  
  //const auto C00 = A(0,0)*B(0,0)+A(0,1)*B(1,0) + A(0,2)*B(2,0);
  //const auto C10 = A(1,0)*B(0,0)+A(1,1)*B(1,0) + A(1,2)*B(2,0);
  //const auto C20 = A(2,0)*B(0,0)+A(2,1)*B(1,0) + A(2,2)*B(2,0);
  //const auto C01 = A(0,0)*B(0,1)+A(0,1)*B(1,1) + A(0,2)*B(2,1);
  //const auto C11 = A(1,0)*B(0,1)+A(1,1)*B(1,1) + A(1,2)*B(2,1);
  //const auto C21 = A(2,0)*B(0,1)+A(2,1)*B(1,1) + A(2,2)*B(2,1);
  //const auto C02 = A(0,0)*B(0,2)+A(0,1)*B(1,2) + A(0,2)*B(2,2);
  //const auto C12 = A(1,0)*B(0,2)+A(1,1)*B(1,2) + A(1,2)*B(2,2);  
  //const auto C22 = A(2,0)*B(0,2)+A(2,1)*B(1,2) + A(2,2)*B(2,2);
  
  return tensor2_9cm< Scal>  {C00, C10, C20, C01, C11, C21, C02, C12, C22 };
}
*/



template <class Scal>
tensor2_9cm< Scal> operator*(Scal alpha, const tensor2_9cm< Scal> &A){
  tensor2_9cm< Scal> B(A);
  return B*=alpha;
}
template <class Scal>
tensor2_9cm< Scal> operator*(const tensor2_9cm< Scal> &A, Scal alpha){
  tensor2_9cm< Scal> B(A);
  return B*=alpha;
}
template <class Scal>
tensor2_9cm< Scal> operator/(Scal alpha, const tensor2_9cm< Scal> &A){
  tensor2_9cm< Scal> B(A);
  return B/=alpha;
}
template <class Scal>
tensor2_9cm< Scal> operator/(const tensor2_9cm< Scal> &A, Scal alpha){
  tensor2_9cm< Scal> B(A);
  return B/=alpha;
}

template <class Scal>
tensor2_9cm< Scal> operator+(const tensor2_9cm< Scal> &A, const tensor2_9cm< Scal> &B){
  tensor2_9cm< Scal> C(A);
  return C+=B;
  /*
  return tensor2_9cm<Scal> {
    A.data[0]+B.data[0],
      A.data[1]+B.data[1],
      A.data[2]+B.data[2],
      A.data[3]+B.data[3],
      A.data[4]+B.data[4],
      A.data[5]+B.data[5],
      A.data[6]+B.data[6],
      A.data[7]+B.data[7],
      A.data[8]+B.data[8]};
  */
  /*tensor2_9cm< Scal> C(A);
  std::transform(B.data.begin(), B.data.end(), C.data.begin(), C.data.begin(), [](const Scal &b, Scal &c){  return c+b; }  );
  return C;
  */
}
template <class Scal>
tensor2_9cm< Scal> operator-(const tensor2_9cm< Scal> &A, const tensor2_9cm< Scal> &B){
  tensor2_9cm< Scal> C(A);
  return C-=B;
}
template <class Scal>
std::ostream & operator << ( std::ostream &out, const tensor2_9cm< Scal> &A  ){
  out << "[ [" << A.data[0] << " " << A.data[3] << " " << A.data[6] << "], ["
      << A.data[1] << " " << A.data[4] << " " << A.data[7] << "], ["
      << A.data[2] << " " << A.data[5] << " " << A.data[8] << "] ]";
  return out;
}



//SYM
template <class Scal>
class tensor2_sym{
  // 0 -> 00
  // 1 -> 11
  // 2 -> 22
  // 3 -> 12=21
  // 4 -> 02 = 20
  // 5 -> 01 = 10

  // 9cm -> sym
  // 0 -> 00    -> 0
  // 1 -> 10=01 -> 5
  // 2 -> 20=02 -> 4
  // 3 -> 01=10 -> 5
  // 4 -> 11    -> 1
  // 5 -> 21=12 -> 3
  // 6 -> 02=20 -> 4
  // 7 -> 12=21 -> 3
  // 8 -> 22    -> 2
 public:
  typedef Scal value_type;
  tensor2_sym() =default;
  tensor2_sym & operator*=(Scal alpha){
    for(auto & a :data) {a*=alpha;}
    return *this;
   }
  tensor2_sym<Scal> & operator+=(const  tensor2_sym<Scal > &A){
    std::transform(A.data.begin(), A.data.end(), data.begin(), data.begin(),[](const Scal &a, const Scal &b){ return a+b; }  );
    return *this;
  }
public :
  std::array< Scal , 6 > data;
};

template <class Scal>
Scal det(const tensor2_sym<Scal> & _A ){
  const auto & A = _A.data;
  return A[0]*( A[1]*A[2] - A[3]*A[3])
    -    A[5]*(A[2]*A[5]-A[3]*A[4])
    +    A[4]*(A[3]*A[5] - A[4]*A[1]);
}


template <class Scal>
tensor2_sym<Scal> invert(const tensor2_sym<Scal> & _A){
  const Scal detA = det(_A);
  const auto & A = _A.data;
  return tensor2_sym<Scal>
    {   ( A[1] * A[2] - A[3] * A[3] )/detA,
	( A[0] * A[2] - A[4] * A[4] )/detA,
	( A[0] * A[1] - A[5] * A[5] )/detA,
	( A[4] * A[5] - A[0] * A[3] )/detA,
	( A[5] * A[3] - A[4] * A[1] )/detA,
	( A[3] * A[4] - A[2] * A[5] )/detA };
}


template <class Scal>
tensor2_9cm<Scal> operator*( const tensor2_9cm<Scal> & _A,  const tensor2_sym<Scal> & _B){
  const auto & A = _A.data;
  const auto & B = _B.data;
  return tensor2_9cm<Scal>{
    A[0] * B[0] + A[3] *B[5] + A[6]*B[4],
      A[1] * B[0] + A[4] *B[5] + A[7]*B[4],
      A[2] * B[0] + A[5] *B[5] + A[8]*B[4],
      A[0] * B[5] + A[3] *B[1] + A[6]*B[3],
      A[1] * B[5] + A[4] *B[1] + A[7]*B[3],
      A[2] * B[5] + A[5] *B[1] + A[8]*B[3],
      A[0] * B[4] + A[3] *B[3] + A[6]*B[2],
      A[1] * B[4] + A[4] *B[3] + A[7]*B[2],
      A[2] * B[4] + A[5] *B[3] + A[8]*B[2] };
}


template <class Scal >
tensor2_sym<Scal > Green_Lagrange_sym ( const   tensor2_9cm<Scal>  &F){
  const Scal E00 = 0.5*(F.data[0]*F.data[0] + F.data[1]*F.data[1] +F.data[2]*F.data[2] -1.);
  const Scal E11 = 0.5*(F.data[3]*F.data[3] + F.data[4]*F.data[4] + F.data[5]*F.data[5] -1.);
  const Scal E22 = 0.5*(F.data[6]*F.data[6] + F.data[7]*F.data[7] + F.data[8]*F.data[8]-1.);
  const Scal E21 = 0.5*(F.data[6]*F.data[3] + F.data[7]*F.data[4] + F.data[8]*F.data[5]);
  const Scal E20 = 0.5*(F.data[6]*F.data[0] + F.data[7]*F.data[1] + F.data[8]*F.data[2]);
  const Scal E10 = 0.5*(F.data[3]*F.data[0] + F.data[4]*F.data[1] + F.data[5]*F.data[2]);
  return tensor2_sym< Scal> {E00, E11, E22, E21, E20, E10 };
}

template <class Scal >
tensor2_sym<Scal > Right_Cauchy_Green_sym ( const   tensor2_9cm<Scal>  &F){
  const Scal C00 = F.data[0]*F.data[0] + F.data[1]*F.data[1] +F.data[2]*F.data[2];
  const Scal C11 = F.data[3]*F.data[3] + F.data[4]*F.data[4] + F.data[5]*F.data[5];
  const Scal C22 = F.data[6]*F.data[6] + F.data[7]*F.data[7] + F.data[8]*F.data[8];
  const Scal C21 = F.data[6]*F.data[3] + F.data[7]*F.data[4] + F.data[8]*F.data[5];
  const Scal C20 = F.data[6]*F.data[0] + F.data[7]*F.data[1] + F.data[8]*F.data[2];
  const Scal C10 = F.data[3]*F.data[0] + F.data[4]*F.data[1] + F.data[5]*F.data[2];
  return tensor2_sym< Scal> {C00, C11, C22, C21, C20, C10 };
}

template <class Scal>
Scal trace( const   tensor2_sym<Scal>  &A ){
  return A.data[0] + A.data[1] + A.data[2];
}

template <class Scal>
// second invariant = trace ATA
Scal I2( const   tensor2_sym<Scal>  &_A ){
  //ATA00 = A00*A00+ A10*A10 + A20*A20
  //ATA11 = A01*A01+ A11*A11 + A21*A21
  //ATA22 = A02*A02+ A12*A12 + A22*A22
  //  A00*A00+A11*A11 +A22*A22 +2*A12*A12+ 2*A20*A20+ 2*A01*A01
  const auto & A = _A.data;
  return A[0]*A[0] + A[1]*A[1]+ A[2]*A[2] +2.*( A[3]*A[3] + A[4]*A[4] + A[5]*A[5]); 
}

template <class Scal>
tensor2_sym< Scal> operator*(Scal alpha, const tensor2_sym< Scal> &A){
  tensor2_sym< Scal> B(A);
  return B*=alpha;
}
template <class Scal>
tensor2_sym< Scal> operator*(const tensor2_sym< Scal> &A, Scal alpha){
  tensor2_sym< Scal> B(A);
  return B*=alpha;
}

template <class Scal>
tensor2_sym< Scal> operator+(const tensor2_sym< Scal> &A, const tensor2_sym< Scal> &B){
  tensor2_sym< Scal> C(A);
  return C+=A;
  /*
  return tensor2_sym<Scal> {
    A.data[0]+B.data[0],
      A.data[1]+B.data[1],
      A.data[2]+B.data[2],
      A.data[3]+B.data[3],
      A.data[4]+B.data[4],
      A.data[5]+B.data[5],
      A.data[6]+B.data[6],
      A.data[7]+B.data[7],
      A.data[8]+B.data[8]};
  */
  /*tensor2_sym< Scal> C(A);
  std::transform(B.data.begin(), B.data.end(), C.data.begin(), C.data.begin(), [](const Scal &b, Scal &c){  return c+b; }  );
  return C;
  */
}

template <class Scal>
std::ostream & operator << ( std::ostream &out, const tensor2_sym< Scal> &A  ){
  out << "[ [" << A.data[0] << " " << A.data[5] << " " << A.data[4] << "], ["
      << A.data[5] << " " << A.data[1] << " " << A.data[3] << "], ["
      << A.data[4] << " " << A.data[3] << " " << A.data[2] << "] ]";
  return out;
}

//ISOTROPE
template <class Scal>
class tensor2_isotrope{
  //00 -> data;
  //01 -> 0
  //02 -> 0
  //10 -> 0
  //11 -> data
  //12 -> 0
  //20 -> 0
  //21 -> 0
  //22 -> data 
 public:
  typedef Scal value_type;
  tensor2_isotrope() =default;
  tensor2_isotrope & operator*=(Scal alpha){
    data*=alpha;
    return *this;
  }
  tensor2_isotrope<Scal> & operator+=(const  tensor2_isotrope<Scal > &A){
    data+=A.data;
    return *this;
  }
public :
  Scal data;
};

template <class Scal>
Scal det(const tensor2_isotrope<Scal> & A ){
  return A.data*A.data*A.data;
}


template <class Scal>
tensor2_isotrope<Scal> invert(const tensor2_isotrope<Scal> & A){
  assert (A.data != 0.);
  return tensor2_isotrope<Scal>{1./A.data};
}


template <class Scal>
tensor2_9cm<Scal> operator*( const tensor2_9cm<Scal> & A,  const tensor2_isotrope<Scal> & B){
  return tensor2_9cm<Scal>{A*B.data};
}

template <class Scal>
tensor2_9cm<Scal> operator*( const tensor2_isotrope<Scal> & A,  const tensor2_9cm<Scal> & B){
  return tensor2_9cm<Scal>{A.data*B};
}

template <class Scal>
tensor2_sym<Scal> operator*( const tensor2_sym<Scal> & A,  const tensor2_isotrope<Scal> & B){
  return tensor2_sym<Scal>{A*B.data};
}

template <class Scal>
tensor2_sym<Scal> operator*( const tensor2_isotrope<Scal> & A,  const tensor2_sym<Scal> & B){
  return tensor2_sym<Scal>{A.data*B};
}

template <class Scal >
tensor2_isotrope<Scal > Green_Lagrange_isotrope ( const   tensor2_isotrope<Scal>  &F){
  return tensor2_isotrope< Scal> {F.data*F.data-1.};
}

template <class Scal >
tensor2_isotrope<Scal > Right_Cauchy_Green_isotrope ( const   tensor2_9cm<Scal>  &F){
  return tensor2_isotrope< Scal> {F.data*F.data };
}

template <class Scal>
Scal trace( const   tensor2_isotrope<Scal>  &A ){
  return 3*A.data;
}


template <class Scal>
tensor2_isotrope< Scal> operator*(Scal alpha, const tensor2_isotrope< Scal> &A){
  tensor2_isotrope< Scal> B(A);
  return B*=alpha;
}
template <class Scal>
tensor2_isotrope< Scal> operator*(const tensor2_isotrope< Scal> &A, Scal alpha){
  tensor2_isotrope< Scal> B(A);
  return B*=alpha;
}

template <class Scal>
tensor2_isotrope< Scal> operator+(const tensor2_isotrope< Scal> &A, const tensor2_isotrope< Scal> &B){
  tensor2_isotrope< Scal> C(A);
  return C+=A;
}

template <class Scal>
tensor2_sym< Scal> operator+(const tensor2_isotrope< Scal> &A, const tensor2_sym< Scal> &B){
  return tensor2_sym< Scal> {B.data[0]+A.data, B.data[1]+A.data, B.data[2]+A.data, B.data[3], B.data[4], B.data[5]   };
}

template <class Scal>
tensor2_9cm< Scal> operator+(const tensor2_isotrope< Scal> &A, const tensor2_9cm< Scal> &B){
  return tensor2_9cm< Scal> {B.data[0]+A.data, B.data[1], B.data[2], B.data[3], B.data[4] + A.data, B.data[5], B.data[6], B.data[7], B.data[8]+A.data   };
}
template <class Scal>
tensor2_9cm< Scal> transpose( const tensor2_9cm < Scal > &_in){
  const auto & in = _in.data;
  return tensor2_9cm< Scal> {
    in[0], in[3], in[6], in[1], in[4], in[7], in[2], in[5], in[8]
    
      };
  
}

template <class Scal>
std::ostream & operator << ( std::ostream &out, const tensor2_isotrope< Scal> &A  ){
  out << "[ [" << A.data << " " << 0      << " " << 0  << "], ["
               << 0      << " " << A.data << " " << 0 << "], ["
               << 0      << " " << 0      << " " << A.data << "] ]";
  return out;
}

template< class Scal>
class mat34cm{
 public:
  std::array<Scal, 12> data;
};

template<class SCAL>
std::ostream & operator << (std::ostream & out, const mat34cm<SCAL > &_A){
  const auto & A = _A.data;
  out << "[ [" << A[0] << " " << A[3] << " " << A[6] << " " << A[9] << "], ["
      << A[1] << " " << A[4] << " " << A[7] << " " << A[10 ] << "], ["
      << A[2] << " " << A[5] << " " << A[8] << " " << A[11] << "] ]";
  return out;
}

template <class Scal>
tensor2_9cm< Scal > apply(const mat34cm<Scal > &_Fop , const mat34cm<Scal > &_x ){
  //Fij = x_ik*Fop_jk
  const auto & x = _x.data;
  const auto & Fop = _Fop.data;
  // std::cout << "x"<< std::endl;
  //std::cout << _x << std::endl;

  //std::cout << "Fop"<< std::endl;
  //std::cout << _Fop << std::endl;
  
  tensor2_9cm<Scal > ret {
    //F00 = x00 Fop00 + X01 Fop01 + X02 Fop02 + X03 Fop03
    x[0] * Fop[0] + x[3]*Fop[3] + x[6]*Fop[6] + x[9]*Fop[9],
    //F10 = x10*Fop00 + x11*Fop01 + x12*Fop02+ x13*Fop03
      x[1] * Fop[0] + x[4]*Fop[3] + x[7]*Fop[6] + x[10]*Fop[9],
      //F20 = x20*Fop00 + x21*Fop01 + x22*Fop02+ x23*Fop03
      x[2] * Fop[0] + x[5]*Fop[3] + x[8]*Fop[6] + x[11]*Fop[9],
      //F01 = x00 Fop10 + X01 Fop11 + X02 Fop12 + X03 Fop13
      x[0] * Fop[1] + x[3]*Fop[4] + x[6]*Fop[7] + x[9]*Fop[10],
      //F11 = x10 Fop10 + X11 Fop11 + X12 Fop12 + X13 Fop13
      x[1] * Fop[1] + x[4]*Fop[4] + x[7]*Fop[7] + x[10]*Fop[10],
      //F21 = x20 Fop10 + X21 Fop11 + X22 Fop12 + X23 Fop13
      x[2] * Fop[1] + x[5]*Fop[4] + x[8]*Fop[7] + x[11]*Fop[10],
      //F02 = x00 Fop20 + X01 Fop21 + X02 Fop22 + X03 Fop23
      x[0] * Fop[2] + x[3]*Fop[5] + x[6]*Fop[8] + x[9]*Fop[11],
      //F12 = x10 Fop20 + X11 Fop21 + X12 Fop22 + X13 Fop23
      x[1] * Fop[2] + x[4]*Fop[5] + x[7]*Fop[8] + x[10]*Fop[11],
      //F22 = x20 Fop20 + X21 Fop21 + X22 Fop22 + X23 Fop23
      x[2] * Fop[2] + x[5]*Fop[5] + x[8]*Fop[8] + x[11]*Fop[11],
      
      
      };
  //std::cout << "F" << std::endl;
  //std::cout << ret << std::endl;
  return ret;

}

template <class Scal>
mat34cm<Scal >  apply(const mat34cm<Scal > &_Fop , const tensor2_9cm<Scal > &_P ){
  //Fij = x_ik*Fop_jk
  const auto & P = _P.data;
  //const auto P = transpose(_P).data;
  const auto & O = _Fop.data;
  // std::cout << "x"<< std::endl;
  //std::cout << _x << std::endl;

  //std::cout << "Fop"<< std::endl;
  //std::cout << _Fop << std::endl;
  
  mat34cm<Scal > ret {
    //Fint_ia = Pij_Oja
    P [0] *O[0] + P[3]*O[1] + P[6]*O[2],
      P [1] *O[0] + P[4]*O[1] + P[7]*O[2],
      P [2] *O[0] + P[5]*O[1] + P[8]*O[2],
      P [0] *O[3] + P[3]*O[4] + P[6]*O[5],
      P [1] *O[3] + P[4]*O[4] + P[7]*O[5],
      P [2] *O[3] + P[5]*O[4] + P[8]*O[5],
      P [0] *O[6] + P[3]*O[7] + P[6]*O[8],
      P [1] *O[6] + P[4]*O[7] + P[7]*O[8],
      P [2] *O[6] + P[5]*O[7] + P[8]*O[8],
      P [0] *O[9] + P[3]*O[10] + P[6]*O[11],
      P [1] *O[9] + P[4]*O[10] + P[7]*O[11],
      P [2] *O[9] + P[5]*O[10] + P[8]*O[11]
      };
  //std::cout << "F" << std::endl;
  //std::cout << ret << std::endl;
  return ret;

}
#endif
