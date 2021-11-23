#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

#include <iostream>
#include <cassert>

#include "espic_type.h"

class ScalarField {
  friend class VectorField;

  public: 
    typedef Real&       reference;
    typedef const Real& const_reference;

    /* Constructors */
    explicit ScalarField(const class Mesh&);

    explicit ScalarField(const class Mesh*);

    /* Copy constructor */
    ScalarField(const ScalarField&);

    /* Destructor */
    ~ScalarField();

    /* Public methods */
    // _ndim = 5 if axi-symmetric
    int ndim()  const { return (_ndim == 5 ? 2 : _ndim); }
    Index dim(const Index) const;
    bool  empty() const { return _size == 0; }
    Index size()  const { return _size; }
    const class Mesh& mesh()  const { return _mesh; }

    /* Assignment operator */
    ScalarField& operator=(const ScalarField&);

    ScalarField& operator=(const Vector&);

    ScalarField& operator=(const Real);

    ScalarField& operator+=(const ScalarField&);

    ScalarField& operator+=(const Vector&);

    ScalarField& operator+=(const Real);

    ScalarField& operator -=(const ScalarField&);

    ScalarField& operator -=(const Vector&);

    ScalarField& operator -=(const Real);

    ScalarField& operator *=(const Real);
    ScalarField& operator /=(const Real);

    ScalarField operator -();     // negative sign, unary operator

    void copy_to(Vector&) const;  // copy data to a Vector

    /* Reference */
          reference operator [] (const Index);
    const_reference operator [] (const Index) const;
          reference operator () (const Index);
    const_reference operator () (const Index) const;
          reference operator () (const Index, const Index);
    const_reference operator () (const Index, const Index) const;
          reference operator () (const Index, const Index, const Index);
    const_reference operator () (const Index, const Index, const Index) const;

    const class Mesh& _mesh;
  private:
    int _ndim;        // 2D or 3D
    Index _dim[3];      // dimension in each direction
    Index _size;        // number of entries
    Real* _data;
};

ScalarField operator+(const ScalarField&, const ScalarField&);
ScalarField operator+(const ScalarField&, const Vector&);
ScalarField operator+(const Vector&, const ScalarField&);
ScalarField operator+(const ScalarField&, const Real);
ScalarField operator+(const Real, const ScalarField&);

ScalarField operator -(const ScalarField&, const ScalarField&);
ScalarField operator -(const ScalarField&, const Vector&);
ScalarField operator -(const Vector&, const ScalarField&);
ScalarField operator -(const ScalarField&, const Real);
ScalarField operator -(const Real, const ScalarField&);

ScalarField operator *(const ScalarField&, const Real);
ScalarField operator *(const Real, const ScalarField&);
ScalarField operator /(const ScalarField&, const Real);

/* inline member functions */
inline Index ScalarField::dim(const Index i) const { 
#ifdef DEBUG
  assert(i < 3);
#endif
  return _dim[i]; 
}

/* Interface methods */
inline ScalarField::reference 
ScalarField::operator[] (const Index n) { 
#ifdef DEBUG
  assert(n < size());
#endif
  return _data[n];
}

inline ScalarField::const_reference 
ScalarField::operator[] (const Index n) const {
#ifdef DEBUG
  assert(n < size());
#endif
  return _data[n];
}

inline ScalarField::reference 
ScalarField::operator () (const Index n) {
#ifdef DEBUG
  assert(n < size());
#endif
  return _data[n];
}

inline ScalarField::const_reference 
ScalarField::operator () (const Index n) const {
#ifdef DEBUG
  assert(n < size());
#endif
  return _data[n];
}

inline ScalarField::reference 
ScalarField::operator () (const Index ix, const Index iy) {
#ifdef DEBUG
  assert(2 == _ndim || 5 == _ndim);
  assert(ix < _dim[0] && iy < _dim[1]);
#endif
  return _data[iy*_dim[0] + ix];
}

inline ScalarField::const_reference 
ScalarField::operator () (const Index ix, const Index iy) const {
#ifdef DEBUG
  assert(2 == _ndim || 5 == _ndim);
  assert(ix < _dim[0] && iy < _dim[1]);
#endif
  return _data[iy*_dim[0] + ix];
}

inline ScalarField::reference 
ScalarField::operator () (const Index ix, const Index iy, const Index iz) {
#ifdef DEBUG
  assert(3 == _ndim);
  assert(ix < _dim[0] && iy < _dim[1] && iz < _dim[2]);
#endif
  return _data[(iz*_dim[1]+iy)*_dim[0]+ix];
}

inline ScalarField::const_reference 
ScalarField::operator () (const Index ix, const Index iy, const Index iz) const {
#ifdef DEBUG
  assert(3 == _ndim);
  assert(ix < _dim[0] && iy < _dim[1] && iz < _dim[2]);
#endif
  return _data[(iz*_dim[1]+iy)*_dim[0]+ix];
}

inline std::ostream& operator<< (std::ostream& os, const ScalarField& f)
{
  for (Index iz = 0; iz < f.dim(2); ++iz) {
    for (Index iy = 0; iy < f.dim(1); ++iy) {
      Index nid = (iz*f.dim(1) + iy)*f.dim(0);

      for (Index ix = 0; ix < f.dim(0)-1; ++ix, ++nid) {
        os << f(nid) << " ";
      }
      os << f(nid) << "\n";
    }
  }

  return os;
}

#endif
