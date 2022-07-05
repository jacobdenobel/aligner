g++  -O3 -Wall  globalAlign.cpp  \
    && ./a.out  \
    &&  tail -n 1 s1.txt | rev  \
    &&  echo ""  \
    &&  tail -n 1 s2.txt | rev \
    &&  echo  ""

