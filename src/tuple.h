#ifndef TUPLE_H
#define TUPLE_H

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

// This is eugene going a little nuts with C++
template <typename TYPE, unsigned SIZE>
class stuple
{
public:
    
    stuple() { }

    stuple(TYPE first, ...)
    {
        data[ 0 ] = first;

        va_list args;
        va_start( args, first );
        
        for( unsigned i = 1; i < SIZE; ++i )
        {
            data[ i ] = va_arg( args, TYPE ); 
        }

        va_end( args );
    }

    stuple(const TYPE array[SIZE]) {
        memcpy(data, array, SIZE * sizeof(TYPE));
    }
    
    stuple(const stuple& other) {
        memcpy(this->data, &other.data, SIZE * sizeof(TYPE));
    }

    stuple& operator=(const stuple& other) {
        if (&other != this)
            memcpy(this->data, &other.data, SIZE * sizeof(TYPE));
        return *this;
    }
    
    stuple& operator=(const TYPE array[SIZE]) {
        memcpy(data, array, SIZE * sizeof(TYPE));
        return *this;
    }

    inline TYPE& operator[](unsigned i) {
        return data[i];
    }

    inline const TYPE& operator[](unsigned i) const {
        return data[i];
    }
private:
    TYPE data[SIZE];
    
};


template <typename TYPE, unsigned SIZE>
std::ostream& operator<<(std::ostream &out, const stuple<TYPE,SIZE> x)
{
    out << '[';
    for (unsigned i=0; i<SIZE; i++)
        out << ' ' << x[i];
    return out << " ]" <<  std::flush;
}

template <typename TYPE, unsigned SIZE>
bool operator == (const stuple<TYPE, SIZE> first, const stuple<TYPE, SIZE> second) {
	return first[0] == second[0] && first[1] == second[1];
}

template <typename TYPE, unsigned SIZE>
bool operator != (const stuple<TYPE, SIZE> first, const stuple<TYPE, SIZE> second) {
	return !(first == second);
}

template <typename TYPE, unsigned SIZE>
bool operator < (const stuple<TYPE, SIZE> first, const stuple<TYPE, SIZE> second) {
	if (first[0] == second[0]) {
		return first[1] < second[1];
	} else {
		return first[0] < second[0];
	}
}

template <typename TYPE, unsigned SIZE>
bool operator <= (const stuple<TYPE, SIZE> first, const stuple<TYPE, SIZE> second) {
	if (first[0] == second[0]) {
		return first[1] <= second[1];
	} else {
		return first[0] < second[0];
	}
}

template <typename TYPE, unsigned SIZE>
bool operator > (const stuple<TYPE, SIZE> first, const stuple<TYPE, SIZE> second) {
	if (first[0] == second[0]) {
		return first[1] > second[1];
	} else {
		return first[0] > second[0];
	}
}

template <typename TYPE, unsigned SIZE>
bool operator >= (const stuple<TYPE, SIZE> first, const stuple<TYPE, SIZE> second) {
	if (first[0] == second[0]) {
		return first[1] >= second[1];
	} else {
		return first[0] > second[0];
	}
}


#endif
