#include <iostream>
#include <cstdio>
#include <limits>
#include <cstdint>
#include <iomanip>

struct TData{
    double key;
    char val[66] = {'\0'};
};

template<class T>
class TVector {
private:
    size_t size_;
    size_t cap_;
    T *data_;
public:
    TVector() :
            size_(0), cap_(0), data_(nullptr) {};

    TVector(double n) :
            size_(n), cap_(n), data_(new T[cap_]) {};

    TVector(double n, T x){
        size_ = n;
        cap_ = n;
        data_ = new T[cap_];
        for (size_t i = 0; i < size_; i++)
            data_[i] = x;
    }

    TVector(const TVector<T>& other){
        if (data_)
            delete[] data_;
        data_ = new T[other.cap_];
        for (size_t i = 0; i < other.size_; ++i) {
            data_[i] = other.data_[i];
        }
        size_ = other.Size();
        cap_ = other.cap_;
    }

    ~TVector(){
        delete[] data_;
    }

    T& operator[] (const int id) const{
        return data_[id];
    }

    TVector<T>& operator= (const TVector<T>& other) {
        if (this != &other) {
            T* tmp = new T[other.size_];
            for (size_t i = 0; i < other.size_; ++i) {
                tmp[i] = other.data_[i];
            }
            delete[] data_;
            data_ = tmp;
            size_ = other.size_;
            cap_ = other.cap_;
        }
        return *this;
    }

    void PushBack(const T& newdata_) {
        if (cap_ == size_){
            cap_ *= 2;
            if (cap_ == 0)
                cap_ = 1;
            T *tmp = new T[cap_];
            for (size_t i = 0; i < size_; ++i) {
                tmp[i] = data_[i];
            }
            delete[] data_;
            data_ = tmp;
        }
        data_[size_++] = newdata_;
    }

    size_t Size() const {
        return size_;
    }
};

/*struct TData{
    double key;
    char val[66] = {'\0'};
}; */

void InsSort(TVector<TData>& vect, double n) {
    for (int i = 1; i < n; i++) {
        TData now = vect[i];
        int j = i - 1;
        while(j >= 0 && vect[j].key > now.key) {
            vect[j + 1] = vect[j];
            --j;
        }
        vect[j + 1] = now;
    }
}

TVector<TData> BucketSort(TVector<TData> const& vect, double n) {
    TVector<TVector<TData> > buckets(n);

    double min_el = std::numeric_limits<double>::max();
    double max_el = std::numeric_limits<double>::min();

    for (int i = 0; i < n; i++) {
        TData elem = vect[i];
        min_el = std::min(min_el, elem.key);
        max_el = std::max(max_el, elem.key);
    }

    double len = (max_el - min_el + 1e-9);

    for (int i = 0; i < n; i++) {
        size_t num = ((vect[i].key - min_el) / len) * n;

        if (num >= (size_t) n)
            num = n - 1;

        buckets[num].PushBack(vect[i]);

    }

    for (int i = 0; i < n; i++) {
        InsSort(buckets[i], buckets[i].Size());
    }

    TVector<TData> res(n);
    int ind = 0;
    for (int i = 0; i < n; i++) {
        for (size_t j = 0; j < buckets[i].Size(); j++) {
            res[ind++] = buckets[i][j];
        }
    }
    return res;
}


    int main() {
        std::ios::sync_with_stdio(false);
        TVector<TData> vect;
        TData in;



        while(std::cin >> in.key >> in.val) {
            vect.PushBack(in);
        }

        TVector<TData> res;
        res = BucketSort(vect, vect.Size());

        for (size_t i = 0; i < res.Size(); i++) {
            std::cout << std::fixed << std::setprecision(6) << res[i].key << '\t' << res[i].val << "\n";
        }
        return 0;
    }
