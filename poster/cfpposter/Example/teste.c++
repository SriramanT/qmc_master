#include <iostream>

template <bool TESTE>
void f (int a)
{
  if(TESTE)
    std::cout << " teste" << a << std::endl;
  if(!TESTE)
    std::cout << " AAA" << a << std::endl;
}


int main() {
  f<true> (1);
  f<false> (2);
  
}

