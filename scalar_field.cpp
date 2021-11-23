#include <cstring>

#include "espic_info.h"
#include "scalar_field.h"
#include "mesh.h"

/* ---------------- Begin Public Methods ---------------- */

/* Constructors */
ScalarField::ScalarField(const Mesh& msh) 
  : _mesh(msh), 
    _ndim(msh.dimension()), 
    _size(msh.num_nodes())
{ 
  _data = new Real [size()];
  memset(_data, 0, size()*sizeof(Real));

  _dim[0] = _mesh.num_nodes(0);
  _dim[1] = _mesh.num_nodes(1);
  _dim[2] = _mesh.num_nodes(2);
}

ScalarField::ScalarField(const Mesh* msh) 
  : ScalarField(*msh)
{
}

/* Copy constructor */
ScalarField::ScalarField(const ScalarField& orig) 
  : _mesh(orig._mesh), 
    _ndim(orig._ndim), 
    _size(orig._size)
{
  _data = new Real [size()];
  memcpy(_data, orig._data, size()*sizeof(Real));

  _dim[0] = orig.dim(0); _dim[1] = orig.dim(1); _dim[2] = orig.dim(2);
}

/* Destructor */
ScalarField::~ScalarField() { 
  delete [] _data;
}


/* Interface */

/* Assignment ScalarField::operator */
ScalarField& ScalarField::operator=(const ScalarField& rhs) {
  if (this == &rhs) return *this;       // standard alias test

  if (&_mesh != &(rhs._mesh)) {
    espic_error("Trying assignment between ScalarFields based on different \"mesh\"");
  }

  memcpy(_data, rhs._data, size()*sizeof(Real));

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator=(const Vector& rhs) {
  if (rhs.size() != size()) {
    espic_error("Trying assignment between ScalarField and Vector with different sizes");
  }

  for (Index n = 0; n < size(); ++n) _data[n] = rhs[n];

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator=(const Real val) {
  for (Index n = 0; n < size(); ++n) _data[n] = val;

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator+=(const ScalarField& rhs) {
  if (rhs.size()!= size()) {
    espic_error("Trying += between two ScalarFields with different sizes");
  }

  for (Index n = 0; n < size(); ++n) _data[n] += rhs._data[n];

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator+=(const Vector& rhs) {
  if (rhs.size() != size()) {
    espic_error("Trying += between ScalarField and Vector with different sizes");
  }

  for (Index n = 0; n < size(); ++n) _data[n] += rhs[n];

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator+=(const Real val) {
  for (Index n = 0; n < size(); ++n) _data[n] += val;

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator -=(const ScalarField& rhs) {
  if (rhs.size() != size()) {
    espic_error("Trying += between two ScalarFields with different sizes");
  }

  for (Index n = 0; n < size(); ++n) _data[n] -= rhs._data[n];

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator -=(const Vector& rhs) {
  if (rhs.size() != size()) {
    espic_error("Trying += between ScalarField and Vector with different sizes");
  }

  for (Index n = 0; n < size(); ++n) _data[n] -= rhs[n];

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator -=(const Real val) {
  for (Index n = 0; n < size(); ++n) _data[n] -= val;

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator *=(const Real val) {
  for (Index n = 0; n < size(); ++n) _data[n] *= val;

  return *this;
}

/* ------------------------------------------------------- */

ScalarField& ScalarField::operator /=(const Real val) {
  for (Index n = 0; n < size(); ++n) _data[n] /= val;

  return *this;
}

/* ------------------------------------------------------- */

ScalarField ScalarField::operator -() { // negative sign, unary operator
  ScalarField ret(*this);
  for (Index n = 0; n < size(); ++n) ret[n] = -ret[n];
  return ret;
}


/* ------------------------------------------------------- */

void ScalarField::copy_to(Vector& v) const {
  if (v.size() != size()) {
    espic_error("Copy ScalarField to Vector with different sizes");
  }

  for (Index n = 0; n < size(); ++n) v[n] = _data[n];
}
/* ----------------- End Public Methods ----------------- */

/* ---------------- Begin Non-Member Functions ---------------- */
ScalarField operator+(const ScalarField& lhs, const ScalarField& rhs)
{
  ScalarField ret(lhs);
  ret += rhs;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator+(const ScalarField& sf, const Vector& vec)
{
  ScalarField ret(sf);
  ret += vec;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator+(const Vector& vec, const ScalarField& sf)
{
  ScalarField ret(sf);
  ret += vec;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator+(const ScalarField& sf, const Real val)
{
  ScalarField ret(sf);
  ret += val;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator+(const Real val, const ScalarField& sf)
{
  ScalarField ret(sf);
  ret += val;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator -(const ScalarField& lhs, const ScalarField& rhs)
{
  ScalarField ret(lhs);
  ret -= rhs;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator -(const ScalarField& sf, const Vector& vec)
{
  ScalarField ret(sf);
  ret -= vec;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator -(const Vector& vec, const ScalarField& sf)
{
  ScalarField ret(sf);    // temporary variable created only once
  ret = vec; ret -= sf;   // no temporary variable created this line
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator -(const ScalarField& sf, const Real val)
{
  ScalarField ret(sf);
  ret -= val;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator -(const Real val, const ScalarField& sf)
{
  ScalarField ret(sf);    // temporary variable created only once
  ret = val; ret -= sf;   // no temporary variable created this line
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator *(const ScalarField& sf, const Real val)
{
  ScalarField ret(sf);    // temporary variable created only once
  ret *= val;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator *(const Real val, const ScalarField& sf)
{
  ScalarField ret(sf);    // temporary variable created only once
  ret *= val;
  return ret;
}

/* ------------------------------------------------------- */

ScalarField operator /(const ScalarField& sf, const Real val)
{
  ScalarField ret(sf);    // temporary variable created only once
  ret /= val;
  return ret;
}
/* ---------------- End Non-Member Functions ---------------- */
