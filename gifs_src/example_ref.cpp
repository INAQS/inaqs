#include <iostream>
#include <string>
#include <vector>


template<typename T>
class VectorReference
{
public:
    VectorReference(T* in_ref) : ref{in_ref} {}
    inline T& vec() {return *ref;}
    inline void update(T* in_ref) {ref = in_ref;}
private:
    T* ref{nullptr};
};


template<typename T>
class UpdatedVectorReference
{
public:
    UpdatedVectorReference(VectorReference<T>& in_ref) : ref{in_ref} {}
    inline T& vec() {return ref.vec();}
private:
    VectorReference<T>& ref{nullptr};
};


int
main()
{
    std::vector<double> va{11231, 2, 3, 4};
    std::vector<double> vb{2, 3, 4};
    VectorReference vref(&va);
    UpdatedVectorReference vvref{vref};

    auto& vec = vvref.vec();

    std::cout << vec[0] << "\n";

    vref.update(&vb);

    vec = vvref.vec();

    std::cout << vec[0] << "\n";

};


