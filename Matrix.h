//
// Created by Dominik on 11.09.26.
//
#ifndef CPP_PROTOTYPE_MATRIX_H
#define CPP_PROTOTYPE_MATRIX_H

#include "Constants.h"

#include <memory>
#include <Accelerate/Accelerate.h>

namespace datastruct {
    template <typename T>
class Matrix {
    public:
        Matrix(T val)
            :
              storage(std::make_unique_for_overwrite<T[]>(elem_cnt)),
              raw_ptr(storage.get()) {
            for (int j = 0; j < ny + 2; j++) {
                for (int i = 0; i < nx + 2; i++) {
                    raw_ptr[j * (nx + 2) + i] = val;
                    if (i == 0) {
                        raw_ptr[j * (nx + 2) + i] = 0.0;
                    }
                    else if (i == nx + 1) {
                        raw_ptr[j * (nx + 2) + i] = 0.0;
                    }
                    if (j == 0) {
                        raw_ptr[j * (nx + 2) + i] = 0.0;
                    }
                    else if (j == ny + 1) {
                        raw_ptr[j * (nx + 2) + i] = 0.0;
                    }
                }
            }
        }
        Matrix()
            :
              storage(std::make_unique_for_overwrite<T[]>(elem_cnt)),
              raw_ptr(storage.get()) {
            for (int j = 0; j < ny + 2; j++) {
                for (int i = 0; i < nx + 2; i++) {
                    raw_ptr[j * (nx + 2) + i] = static_cast<T>(0.0);
                }
            }
        }

        Matrix(T* external_ptr)
            :
              storage(nullptr),
              raw_ptr(external_ptr) {}

        // Move semantics (Rule of 5)
        Matrix(Matrix&&) noexcept = default;
        Matrix& operator=(Matrix&&) noexcept = default;
        Matrix(const Matrix&) = delete;
        Matrix& operator=(const Matrix&) = delete;

        void clone(const Matrix &m) {
            cblas_dcopy(elem_cnt, m.data(), 1, this->data(), 1);
        }

        // 2D access: row j, col i
        T& operator()(const size_t j, const size_t i) noexcept {
            return raw_ptr[j * (nx + 2) + i];
        }

        const T& operator()(const size_t j, const size_t i) const noexcept {
            return raw_ptr[j * (nx + 2) + i];
        }

        void subtract(const Matrix &m, double scalar_a, double scalar_b) {
            catlas_daxpby(elem_cnt,  -scalar_a, m.data(), 1, scalar_b, this->data(), 1);
        }

        void subtract(const Matrix &m, double scalar_a) {
            catlas_daxpby(elem_cnt,  -scalar_a, m.data(), 1, 1.0, this->data(), 1);
        }

        void subtract(const Matrix &m) {
            catlas_daxpby(elem_cnt,  -1.0, m.data(), 1, 1.0, this->data(), 1);
        }

        int maxIdx() {
            return cblas_idamax(elem_cnt, this->data(), 1);
        }

        int minIdx() {
            return 0;
        }

        double dot(const Matrix &m) {
            return cblas_ddot(elem_cnt, this->data(), 1, m.data(), 1);
        }

        std::string printMatrix() {
            std::string outstr = "[";
            for (int j= 1; j < ny + 1; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    outstr.append(std::format("{:.{}g}", raw_ptr[j * (nx + 2) + i], std::numeric_limits<double>::max_digits10));
                    outstr.append(", ");
                }
            }
            outstr.replace(outstr.size() - 2, 2, "]");
            return outstr;
        }

        double norm() {
            return cblas_dnrm2(elem_cnt, this->data(), 1);
        }

        double mean() {
            return cblas_ddot(elem_cnt, this->data(), 1, ones().data(), 1)/static_cast<double>(nx*ny);
        }

        void add(const Matrix &m) {
            catlas_daxpby(elem_cnt, 1.0, m.data(), 1, 1.0, this->data(), 1);
        }

        void add(const Matrix &m, double scalar_a) {
            catlas_daxpby(elem_cnt, scalar_a, m.data(), 1, 1.0, this->data(), 1);
        }

        void add(const Matrix &m, double scalar_a, double scalar_b) {
            catlas_daxpby(elem_cnt, scalar_a, m.data(), 1, scalar_b, this->data(), 1);
        }

        static const Matrix& ones() {
            static const Matrix instance(1.0);
            return instance;
        }

        __attribute__((noinline, used))
        std::vector<std::vector<double>> collectAsVector() const {
            std::vector<std::vector<double>> res;
            for (int j = 0; j < ny + 2; j++) {
                std::vector<double> help(nx + 2);
                for (int i = 0; i < nx + 2; i++) {
                    help[i] = (*this)(j, i);
                }
                res.push_back(help);
            }
            return res;
        }

        // 1D flat access
        T& operator()(size_t i) noexcept { return raw_ptr[i]; }
        const T& operator()(size_t i) const noexcept { return raw_ptr[i]; }

        // Raw pointers for vectorization/SIMD
        T* data() noexcept { return raw_ptr; }
        const T* data() const noexcept { return raw_ptr; }

    private:
        std::unique_ptr<T[]> storage;
        T* raw_ptr;
};
}
#endif //CPP_PROTOTYPE_MATRIX_H