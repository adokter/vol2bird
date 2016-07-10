/* --------------------------------------------------------------------
Copyright (C) 2011 Swedish Meteorological and Hydrological Institute, SMHI,

This file is part of beamb.

beamb is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

beamb is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with beamb.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------*/
/**
 * Python version of the beam blockage
 * @file
 * @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
 * @date 2011-11-14
 */
#ifndef PYVOL2BIRD_H
#define PYVOL2BIRD_H
#include "libvol2bird.h"

/**
 * The beam blockage
 */
typedef struct {
   PyObject_HEAD /*Always have to be on top*/
   vol2bird_t* v2b;  /**< the native object */
} PyVol2Bird;

#define PyVol2Bird_Type_NUM 0                              /**< index of type */

#define PyVol2Bird_GetNative_NUM 1                         /**< index of GetNative*/
#define PyVol2Bird_GetNative_RETURN vol2bird_t*         /**< return type for GetNative */
#define PyVol2Bird_GetNative_PROTO (PyVol2Bird*)        /**< arguments for GetNative */

#define PyVol2Bird_New_NUM 2                               /**< index of New */
#define PyVol2Bird_New_RETURN PyVol2Bird*              /**< return type for New */
//#define PyVol2Bird_New_PROTO (vol2bird_t*) /**< arguments for New */
#define PyVol2Bird_New_PROTO (PolarVolume_t* volume) /**< arguments for New */

#define PyVol2Bird_API_pointers 3                          /**< number of type and function pointers */

#ifdef PYVOL2BIRD_MODULE
/** Forward declaration of type */
extern PyTypeObject PyVol2Bird_Type;

/** Checks if the object is a PyVol2Bird or not */
#define PyVol2Bird_Check(op) ((op)->ob_type == &PyVol2Bird_Type)

/** Forward declaration of PyVol2Bird_GetNative */
static PyVol2Bird_GetNative_RETURN PyVol2Bird_GetNative PyVol2Bird_GetNative_PROTO;

/** Forward declaration of PyVol2Bird_New */
static PyVol2Bird_New_RETURN PyVol2Bird_New PyVol2Bird_New_PROTO;

#else
/** Pointers to types and functions */
static void **PyVol2Bird_API;

/**
 * Returns a pointer to the internal vol2bird, remember to release the reference
 * when done with the object. (RAVE_OBJECT_RELEASE).
 */
#define PyVol2Bird_GetNative \
  (*(PyVol2Bird_GetNative_RETURN (*)PyVol2Bird_GetNative_PROTO) PyVol2Bird_API[PyVol2Bird_GetNative_NUM])

/**
 * Creates a new vol2bird instance. Release this object with Py_DECREF. If a Vol2Bird_t instance is
 * provided and this instance already is bound to a python instance, this instance will be increfed and
 * returned.
 * @param[in] beamb - the Vol2Bird_t instance.
 * @returns the PyVol2Bird instance.
 */
#define PyVol2Bird_New \
  (*(PyVol2Bird_New_RETURN (*)PyVol2Bird_New_PROTO) PyVol2Bird_API[PyVol2Bird_New_NUM])

/**
 * Checks if the object is a python beam blockage instance
 */
#define PyVol2Bird_Check(op) \
   ((op)->ob_type == (PyTypeObject *)PyVol2Bird_API[PyVol2Bird_Type_NUM])

/**
 * Imports the PyVol2Bird module (like import _pyvol2bird in python).
 */
static int
import_pyvol2bird(void)
{
  PyObject *module;
  PyObject *c_api_object;

  module = PyImport_ImportModule("_pyvol2bird");
  if (module == NULL) {
    return -1;
  }

  c_api_object = PyObject_GetAttrString(module, "_C_API");
  if (c_api_object == NULL) {
    Py_DECREF(module);
    return -1;
  }
  if (PyCObject_Check(c_api_object)) {
    PyVol2Bird_API = (void **)PyCObject_AsVoidPtr(c_api_object);
  }
  Py_DECREF(c_api_object);
  Py_DECREF(module);
  return 0;
}

#endif

#endif /* PYVOL2BIRD_H */
