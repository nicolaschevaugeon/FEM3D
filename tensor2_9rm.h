template <class Scal >
tensor2_9rm<Scal > Green_Lagrange ( const   tensor2_9rm<Scal>  &F){
  const auto E00 = 0.5*(F(0,0)*F(0,0) + F(1,0)*F(1,0)+ F(2,0)*F(2,0) -1.);
  const auto E01 = 0.5*(F(0,0)*F(0,1) + F(1,0)*F(1,1)+ F(2,0)*F(2,1) -1.);
  const auto E02 = 0.5*(F(0,0)*F(0,2) + F(1,0)*F(1,2)+ F(2,0)*F(2,2) -1.);
  const auto E10 = 0.5*(F(0,1)*F(0,0) + F(1,1)*F(1,0)+ F(2,1)*F(2,0) -1.);
  const auto E11 = 0.5*(F(0,1)*F(0,1) + F(1,1)*F(1,1)+ F(2,1)*F(2,1) -1.);
  const auto E12 = 0.5*(F(0,1)*F(0,2) + F(1,1)*F(1,2)+ F(2,1)*F(2,2) -1.);
  const auto E20 = 0.5*(F(0,2)*F(0,0) + F(1,2)*F(1,0)+ F(2,2)*F(2,0) -1.);
  const auto E21 = 0.5*(F(0,2)*F(0,1) + F(1,2)*F(1,1)+ F(2,2)*F(2,1) -1.);
  const auto E22 = 0.5*(F(0,2)*F(0,2) + F(1,2)*F(1,2)+ F(2,2)*F(2,2) -1.);
  //  return tensor2_9rm< double> (std::array<Scal, 9> {E00, E01, E02, E10, E11, E12, E20, E21, E22 });
  return tensor2_9rm< double>{E00, E01, E02, E10, E11, E12, E20, E21, E22 };

}


template <class Scal>
class tensor2_9rm{
  // template < class T> tensor2_9rm( const T &in):data(std::forward(in)){}
public:
  tensor2_9rm() = default;
  tensor2_9rm(const tensor2_9rm<Scal> & ) = default;
  
  /* tensor2_9rm( std::initializer_list<Scal > _data){
    //static_assert(_data.size()< 9, "gg");
    std::copy(_data.begin(), _data.end(), data.begin());
    };*/
  
  tensor2_9rm(tensor2_9rm<Scal> && ) = default;
  tensor2_9rm & operator = (tensor2_9rm<Scal> && ) = default;
  //tensor2_9rm( std::array< Scal, 9> &&_data):data(std::move(_data)){}
  Scal &operator()(size_t i, size_t j ) {return data[3*i+j];}
  const Scal &operator()(size_t i, size_t j ) const {return data[3*i+j];}
  tensor2_9rm & operator*=(Scal alpha){
    for(auto & a :data) {a*=alpha;}
    return *this;
  }
  tensor2_9rm<Scal> & operator+=(const  tensor2_9rm<Scal > &A){
    std::transform(A.data.begin(), A.data.end(), data.begin(), data.begin(),[](const Scal &a, const Scal &b){ return a+b; }  );
    return *this;
  }
public :
  std::array< Scal , 9 > data;
  // std::vector<double > datareturn data[i+j];;
};


template <class Scal >
tensor2_9rm<Scal > operator* ( const   tensor2_9rm<Scal>  &A, const   tensor2_9rm<Scal>  &B){
  /*const auto C00= A.data[0]*B.data[0] + A.data[1]*B.data[3] + A.data[2]*B.data[6]; 
  const auto C01= A.data[0]*B.data[1] + A.data[1]*B.data[4] + A.data[2]*B.data[7];
  const auto C02= A.data[0]*B.data[2] + A.data[1]*B.data[5] + A.data[2]*B.data[8];
  const auto C10= A.data[3]*B.data[0] + A.data[4]*B.data[3] + A.data[5]*B.data[6]; 
  const auto C11= A.data[3]*B.data[1] + A.data[4]*B.data[4] + A.data[5]*B.data[7];
  const auto C12= A.data[3]*B.data[2] + A.data[4]*B.data[5] + A.data[5]*B.data[8];
  const auto C20= A.data[6]*B.data[0] + A.data[7]*B.data[3] + A.data[8]*B.data[6]; 
  const auto C21= A.data[6]*B.data[1] + A.data[7]*B.data[4] + A.data[8]*B.data[7];
  const auto C22= A.data[6]*B.data[2] + A.data[7]*B.data[5] + A.data[8]*B.data[8];*/
  
  const auto C00 = A(0,0)*B(0,0)+A(0,1)*B(1,0) + A(0,2)*B(2,0);
  const auto C01 = A(0,0)*B(0,1)+A(0,1)*B(1,1) + A(0,2)*B(2,1);
  const auto C02 = A(0,0)*B(0,2)+A(0,1)*B(1,2) + A(0,2)*B(2,2);
  const auto C10 = A(1,0)*B(0,0)+A(1,1)*B(1,0) + A(1,2)*B(2,0);
  const auto C11 = A(1,0)*B(0,1)+A(1,1)*B(1,1) + A(1,2)*B(2,1);
  const auto C12 = A(1,0)*B(0,2)+A(1,1)*B(1,2) + A(1,2)*B(2,2);
  const auto C20 = A(2,0)*B(0,0)+A(2,1)*B(1,0) + A(2,2)*B(2,0);
  const auto C21 = A(2,0)*B(0,1)+A(2,1)*B(1,1) + A(2,2)*B(2,1);
  const auto C22 = A(2,0)*B(0,2)+A(2,1)*B(1,2) + A(2,2)*B(2,2);
  // return tensor2_9rm< Scal> (std::array<Scal, 9> {C00, C01, C02, C10, C11, C12, C20, C21, C22 });

  return tensor2_9rm< Scal> {C00, C01, C02, C10, C11, C12, C20, C21, C22 };
}



template <class Scal>
tensor2_9rm< Scal> operator*(Scal alpha, const tensor2_9rm< Scal> &A){
  tensor2_9rm< Scal> B(A);
  return B*=alpha;
}
template <class Scal >
tensor2_9rm< Scal> operator*(const tensor2_9rm< Scal> &A, Scal alpha){
  tensor2_9rm< Scal> B(A);
  return B*=alpha;
}

template <class Scal>
tensor2_9rm< Scal> operator+(const tensor2_9rm< Scal> &A, const tensor2_9rm< Scal> &B){
  tensor2_9rm< Scal> C(A);
  return C+=A;
}


template <class Scal>
Scal trace( const   tensor2_9rm<Scal>  &A ){
  return A(0,0)+ A(1,1)+A(2,2);
}
