#pragma once

#include <cassert>

template <typename T>
class genlist {
private:
    int Size;
    T* data;

public:
    const int& size;

    genlist()
        : Size(0), data(nullptr), size(Size)
    {
    }

    genlist(const genlist& other)
        : Size(other.Size), data(new T[other.Size]()), size(Size)
    {
        for (int i = 0; i < Size; i++) {
            data[i] = other.data[i];
        }
    }

    genlist(genlist&& other) = delete;

    ~genlist()
    {
        delete[] data;
    }

    genlist& operator=(const genlist& other)
    {
        if (this == &other) {
            return *this;
        }

        T* newdata = new T[other.Size]();

        for (int i = 0; i < other.Size; i++) {
            newdata[i] = other.data[i];
        }

        delete[] data;

        Size = other.Size;
        data = newdata;

        return *this;
    }

    genlist& operator=(genlist&& other) = delete;

    T& operator[](int i)
    {
        assert(i >= 0 && i < Size);
        return data[i];
    }

    const T& operator[](int i) const
    {
        assert(i >= 0 && i < Size);
        return data[i];
    }

    void add(const T& item)
    {
        T* newdata = new T[Size + 1]();

        for (int i = 0; i < Size; i++) {
            newdata[i] = data[i];
        }

        newdata[Size] = item;

        delete[] data;

        data = newdata;
        Size++;
    }

    void remove(int index)
    {
        assert(index >= 0 && index < Size);

        T* newdata = new T[Size - 1]();

        for (int old_i = 0, new_i = 0; old_i < Size; old_i++) {
            if (old_i != index) {
                newdata[new_i] = data[old_i];
                new_i++;
            }
        }

        delete[] data;

        data = newdata;
        Size--;
    }
};