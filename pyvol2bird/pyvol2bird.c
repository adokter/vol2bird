/* --------------------------------------------------------------------
Copyright (C) 2016 Swedish Meteorological and Hydrological Institute, SMHI,

This file is part of vol2bird.

vol2bird is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

vol2bird is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with vol2bird.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------*/
/**
 * Python API to the main vol2bird function
 * @file
 * @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
 * @author Adriaan Dokter (University of Amsterdam, UvA)
 * @author Jurriaan Spaaks, Lourens Veen (Netherlands eScience centre, NLeSC)
 * @date 2016-06-14
 */
#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <string.h>


#define PYVOL2BIRD_MODULE   /**< to get correct part in pyvol2bird */
#include "pyvol2bird.h"

#include "pyverticalprofile.h"
#include "pypolarvolume.h"
#include "pyrave_debug.h"
#include "rave_alloc.h"

/**
 * Debug this module
 */
PYRAVE_DEBUG_MODULE("_pyvol2bird");

/**
 * Sets a python exception and goto tag
 */
#define raiseException_gotoTag(tag, type, msg) \
{PyErr_SetString(type, msg); goto tag;}

/**
 * Sets python exception and returns NULL
 */
#define raiseException_returnNULL(type, msg) \
{PyErr_SetString(type, msg); return NULL;}

/**
 * Error object for reporting errors to the python interpreeter
 */
static PyObject *ErrorObject;

/// --------------------------------------------------------------------
/// Vol2Bird
/// --------------------------------------------------------------------
/*@{ Vol2Bird */
/**
 * Returns the native Vol2Bird_t instance.
 * @param[in] vol - the python Vol2Bird instance
 * @returns the native Vol2Bird_t instance.
 */
static vol2bird_t*
PyVol2Bird_GetNative(PyVol2Bird* v2b)
{
  RAVE_ASSERT((v2b != NULL), "v2b == NULL");
  return RAVE_OBJECT_COPY(v2b->v2b);
}

static PyVol2Bird* PyVol2Bird_New(PolarVolume_t* volume)
{
  PyVol2Bird* result = NULL;
  vol2bird_t* alldata = NULL;

  alldata = malloc(sizeof(vol2bird_t));
  if (alldata == NULL) {
    RAVE_CRITICAL0("Failed to allocate memory for Vol2Bird.");
    raiseException_returnNULL(PyExc_MemoryError, "Failed to allocate memory for Vol2Bird.");
  }

  if (vol2birdLoadConfig(alldata)) {
    raiseException_returnNULL(PyExc_ValueError, "vol2birdLoadConfig did not complete successfully.");
  }

  int initSuccessful = vol2birdSetUp(volume, alldata) == 0;
  if (initSuccessful == FALSE) {
     raiseException_returnNULL(PyExc_ValueError, "vol2birdSetUp did not complete successfully.");
  }
  RAVE_OBJECT_RELEASE(volume);
  result = PyObject_NEW(PyVol2Bird, &PyVol2Bird_Type);
  if (result != NULL) {
    result->v2b = alldata;
  }
  return result;
}

/**
 * Deallocates the beam blockage
 * @param[in] obj the object to deallocate.
 */
static void _pyvol2bird_dealloc(PyVol2Bird* obj)
{
  /*Nothing yet*/
  if (obj == NULL) {
    return;
  }
  //PYRAVE_DEBUG_OBJECT_DESTROYED;

// tear down vol2bird, give memory back
// fairly easy to remove cfg from vol2birdTearDown, to have a single argument destrocutor
// suggested replacement: vol2birdTearDown(cfg, obj->v2b);
  vol2birdTearDown(obj->v2b);
  free(obj->v2b);
  PyObject_Del(obj);
}

/**
 * Creates a new instance of the vol2bird
 * @param[in] self this instance.
 * @param[in] args arguments for creation or a beam blockage
 * @return the object on success, otherwise NULL
 */
static PyObject* _pyvol2bird_new(PyObject* self, PyObject* args)
{
  PyObject* pyin = NULL;
  if (!PyArg_ParseTuple(args, "O", &pyin)) {
    return NULL;
  }

  if (!PyPolarVolume_Check(pyin)) {
    raiseException_returnNULL(PyExc_ValueError, "First argument should be a Polar Scan");
  }

  PolarVolume_t* pvol = PyPolarVolume_GetNative((PyPolarVolume*)pyin);
  return (PyObject*)PyVol2Bird_New(pvol);
}

/**
 * Returns the blockage for the provided scan given gaussian limit.
 * @param[in] self - self
 * @param[in] args - the arguments (PyPolarVolume, double (Limit of Gaussian approximation of main lobe))
 * @return the PyRaveField on success otherwise NULL
 */
static PyObject* _pyvol2bird_vol2bird(PyVol2Bird* self, PyObject* args)
{
  PyObject* pyin = NULL;
  PyObject* result = NULL;

  if (!PyArg_ParseTuple(args, "O", &pyin)) {
    return NULL;
  }

  if (!PyPolarVolume_Check(pyin)) {
    raiseException_returnNULL(PyExc_ValueError, "First argument should be a Polar Scan");
  }

  vol2birdCalcProfiles(self->v2b);
  mapDataToRave(((PyPolarVolume*)pyin)->pvol, self->v2b);
  // do we need a copy of vp?
  if (self->v2b->vp != NULL) {
    result = (PyObject*)PyVerticalProfile_New(self->v2b->vp);
  }
  return result;
}

/**
 * All methods a ropo generator can have
 */
static struct PyMethodDef _pyvol2bird_methods[] =
{
  {"misc_vol2birdSuccessful", NULL},
  {"constants_cellDbzMin", NULL},
  {"vol2bird", (PyCFunction)_pyvol2bird_vol2bird, 1},
  {NULL, NULL} /* sentinel */
};

/**
 * Returns the specified attribute in the beam blockage
 */
static PyObject* _pyvol2bird_getattr(PyVol2Bird* self, char* name)
{
  PyObject* res = NULL;

  if (strcmp("misc_vol2birdSuccessful", name) == 0) {
    return PyInt_FromLong(self->v2b->misc.vol2birdSuccessful);
  } else if(strcmp("constants_cellDbzMin", name) == 0) {
    return PyFloat_FromDouble(self->v2b->constants.cellDbzMin);
  }

  res = Py_FindMethod(_pyvol2bird_methods, (PyObject*) self, name);
  if (res)
    return res;

  PyErr_Clear();
  PyErr_SetString(PyExc_AttributeError, name);
  return NULL;
}

/**
 * Returns the specified attribute in the polar volume
 */
static int _pyvol2bird_setattr(PyVol2Bird* self, char* name, PyObject* val)
{
  int result = -1;
  if (name == NULL) {
    goto done;
  }

  if (strcmp("misc_vol2birdSuccessful", name) == 0) {
    if (PyFloat_Check(val)) {
      self->v2b->misc.vol2birdSuccessful = (int)PyFloat_AsDouble(val);
    } else if (PyLong_Check(val)) {
        self->v2b->misc.vol2birdSuccessful = (int)PyLong_AsDouble(val);
    } else if (PyInt_Check(val)) {
        self->v2b->misc.vol2birdSuccessful = (int)PyInt_AsLong(val);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "misc_vol2birdSuccessful must be number")
    }
  } else if (strcmp("constants_cellDbzMin", name) == 0) {
    if (PyFloat_Check(val)) {
      self->v2b->constants.cellDbzMin = PyFloat_AsDouble(val);
    } else if (PyLong_Check(val)) {
        self->v2b->constants.cellDbzMin = PyLong_AsDouble(val);
    } else if (PyInt_Check(val)) {
        self->v2b->constants.cellDbzMin = (double)PyInt_AsLong(val);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "constants_cellDbzMin must be number")
    }
  } else {
    raiseException_gotoTag(done, PyExc_AttributeError, name);
  }

  result = 0;
done:
  return result;
}

/*@} End of Fmi Image */

/// --------------------------------------------------------------------
/// Type definitions
/// --------------------------------------------------------------------
/*@{ Type definitions */
PyTypeObject PyVol2Bird_Type =
{
  PyObject_HEAD_INIT(NULL)0, /*ob_size*/
  "Vol2BirdCore", /*tp_name*/
  sizeof(PyVol2Bird), /*tp_size*/
  0, /*tp_itemsize*/
  /* methods */
  (destructor)_pyvol2bird_dealloc, /*tp_dealloc*/
  0, /*tp_print*/
  (getattrfunc)_pyvol2bird_getattr, /*tp_getattr*/
  (setattrfunc)_pyvol2bird_setattr, /*tp_setattr*/
  0, /*tp_compare*/
  0, /*tp_repr*/
  0, /*tp_as_number */
  0,
  0, /*tp_as_mapping */
  0 /*tp_hash*/
};
/*@} End of Type definitions */

/*@{ Functions */

/*@} End of Functions */

/*@{ Module setup */
static PyMethodDef functions[] = {
  {"new", (PyCFunction)_pyvol2bird_new, 1},
  {NULL,NULL} /*Sentinel*/
};

PyMODINIT_FUNC
init_pyvol2bird(void)
{
  PyObject *module=NULL,*dictionary=NULL;
  static void *PyVol2Bird_API[PyVol2Bird_API_pointers];
  PyObject *c_api_object = NULL;
  PyVol2Bird_Type.ob_type = &PyType_Type;

  module = Py_InitModule("_pyvol2bird", functions);
  if (module == NULL) {
    return;
  }
  PyVol2Bird_API[PyVol2Bird_Type_NUM] = (void*)&PyVol2Bird_Type;
  PyVol2Bird_API[PyVol2Bird_GetNative_NUM] = (void *)PyVol2Bird_GetNative;
  PyVol2Bird_API[PyVol2Bird_New_NUM] = (void*)PyVol2Bird_New;

  c_api_object = PyCObject_FromVoidPtr((void *)PyVol2Bird_API, NULL);

  if (c_api_object != NULL) {
    PyModule_AddObject(module, "_C_API", c_api_object);
  }

  dictionary = PyModule_GetDict(module);
  ErrorObject = PyString_FromString("_pyvol2bird.error");

  if (ErrorObject == NULL || PyDict_SetItemString(dictionary, "error", ErrorObject) != 0) {
    Py_FatalError("Can't define _pyvol2bird.error");
  }
  import_pypolarvolume();
  import_pyverticalprofile();
  PYRAVE_DEBUG_INITIALIZE;
}
/*@} End of Module setup */
