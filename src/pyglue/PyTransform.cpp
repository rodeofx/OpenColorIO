/*
Copyright (c) 2003-2010 Sony Pictures Imageworks Inc., et al.
All Rights Reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of Sony Pictures Imageworks nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <OpenColorIO/OpenColorIO.h>

#include "PyTransform.h"
#include "PyUtil.h"

#include <sstream>

OCIO_NAMESPACE_ENTER
{
    namespace
    {
        PyOCIO_Transform *  PyTransform_New(ConstTransformRcPtr transform)
        {
            if (!transform)
            {
                return 0x0;
            }
            
            PyOCIO_Transform * pyobj = 0x0;
            
            if(ConstAllocationTransformRcPtr allocationTransform = \
                DynamicPtrCast<const AllocationTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_AllocationTransformType);
            }
            else if(ConstCDLTransformRcPtr cdlTransform = \
                DynamicPtrCast<const CDLTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_CDLTransformType);
            }
            else if(ConstColorSpaceTransformRcPtr colorSpaceTransform = \
                DynamicPtrCast<const ColorSpaceTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_ColorSpaceTransformType);
            }
            else if(ConstDisplayTransformRcPtr displayTransform = \
                DynamicPtrCast<const DisplayTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_DisplayTransformType);
            }
            else if(ConstExponentTransformRcPtr exponentTransform = \
                DynamicPtrCast<const ExponentTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_ExponentTransformType);
            }
            else if(ConstFileTransformRcPtr fileTransform = \
                DynamicPtrCast<const FileTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_FileTransformType);
            }
            else if(ConstGroupTransformRcPtr groupTransform = \
                DynamicPtrCast<const GroupTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_GroupTransformType);
            }
            else if(ConstMatrixTransformRcPtr matrixTransform = \
                DynamicPtrCast<const MatrixTransform>(transform))
            {
                pyobj = PyObject_New(PyOCIO_Transform,
                    (PyTypeObject * ) &PyOCIO_MatrixTransformType);
            }
            
            return pyobj;
        }
    }
    
    PyObject * BuildConstPyTransform(ConstTransformRcPtr transform)
    {
        if (!transform)
        {
            Py_RETURN_NONE;
        }
        
        PyOCIO_Transform * pyobj = PyTransform_New(transform);
        
        if(!pyobj)
        {
            std::ostringstream os;
            os << "Unknown transform type for BuildConstPyTransform.";
            throw Exception(os.str().c_str());
        }
        
        pyobj->constcppobj = new ConstTransformRcPtr();
        pyobj->cppobj = new TransformRcPtr();
        
        *pyobj->constcppobj = transform;
        pyobj->isconst = true;
        
        return (PyObject *) pyobj;
    }
    
    PyObject * BuildEditablePyTransform(TransformRcPtr transform)
    {
        if (!transform)
        {
            Py_RETURN_NONE;
        }
        
        PyOCIO_Transform * pyobj = PyTransform_New(transform);
        
        pyobj->constcppobj = new ConstTransformRcPtr();
        pyobj->cppobj = new TransformRcPtr();
        
        *pyobj->cppobj = transform;
        pyobj->isconst = false;
        
        return (PyObject *) pyobj;
    }
    
    bool IsPyTransform(PyObject * pyobject)
    {
        if(!pyobject) return false;
        return PyObject_TypeCheck(pyobject, &PyOCIO_TransformType);
    }
    
    bool IsPyTransformEditable(PyObject * pyobject)
    {
        if(!IsPyTransform(pyobject)) return false;
        
        PyOCIO_Transform * pyobj = reinterpret_cast<PyOCIO_Transform *> (pyobject);
        return (!pyobj->isconst);
    }
    
    TransformRcPtr GetEditableTransform(PyObject * pyobject)
    {
        if(!IsPyTransform(pyobject))
        {
            throw Exception("PyObject must be an OCIO.Transform.");
        }
        
        PyOCIO_Transform * pytransform = reinterpret_cast<PyOCIO_Transform *> (pyobject);
        if(!pytransform->isconst && pytransform->cppobj)
        {
            return *pytransform->cppobj;
        }
        
        throw Exception("PyObject must be an editable OCIO.Transform.");
    }
    
    ConstTransformRcPtr GetConstTransform(PyObject * pyobject, bool allowCast)
    {
        if(!IsPyTransform(pyobject))
        {
            throw Exception("PyObject must be an OCIO.Transform.");
        }
        
        PyOCIO_Transform * pytransform = reinterpret_cast<PyOCIO_Transform *> (pyobject);
        if(pytransform->isconst && pytransform->constcppobj)
        {
            return *pytransform->constcppobj;
        }
        
        if(allowCast && !pytransform->isconst && pytransform->cppobj)
        {
            return *pytransform->cppobj;
        }
        
        throw Exception("PyObject must be a valid OCIO.Transform.");
    }
    
    
    
    ///////////////////////////////////////////////////////////////////////////
    
    
    bool AddTransformObjectToModule( PyObject* m )
    {
        PyOCIO_TransformType.tp_new = PyType_GenericNew;
        if ( PyType_Ready(&PyOCIO_TransformType) < 0 ) return false;
        
        Py_INCREF( &PyOCIO_TransformType );
        PyModule_AddObject(m, "Transform",
                (PyObject *)&PyOCIO_TransformType);
        
        return true;
    }
    
    
    namespace
    {
        int PyOCIO_Transform_init( PyOCIO_Transform * self, PyObject * args, PyObject * kwds );
        void PyOCIO_Transform_delete( PyOCIO_Transform * self, PyObject * args );
        
        PyObject * PyOCIO_Transform_isEditable( PyObject * self );
        PyObject * PyOCIO_Transform_createEditableCopy( PyObject * self );
        PyObject * PyOCIO_Transform_getDirection( PyObject * self );
        PyObject * PyOCIO_Transform_setDirection( PyObject * self,  PyObject *args );
        
        ///////////////////////////////////////////////////////////////////////
        ///
        
        PyMethodDef PyOCIO_Transform_methods[] = {
            {"isEditable", (PyCFunction) PyOCIO_Transform_isEditable, METH_NOARGS, "" },
            {"createEditableCopy", (PyCFunction) PyOCIO_Transform_createEditableCopy, METH_NOARGS, "" },
            {"getDirection", (PyCFunction) PyOCIO_Transform_getDirection, METH_NOARGS, "" },
            {"setDirection", PyOCIO_Transform_setDirection, METH_VARARGS, "" },
            
            {NULL, NULL, 0, NULL}
        };
    }
    
    ///////////////////////////////////////////////////////////////////////////
    ///
    
    PyTypeObject PyOCIO_TransformType = {
        PyObject_HEAD_INIT(NULL)
        0,                                          //ob_size
        "OCIO.Transform",                           //tp_name
        sizeof(PyOCIO_Transform),                   //tp_basicsize
        0,                                          //tp_itemsize
        (destructor) PyOCIO_Transform_delete,        //tp_dealloc
        0,                                          //tp_print
        0,                                          //tp_getattr
        0,                                          //tp_setattr
        0,                                          //tp_compare
        0,                                          //tp_repr
        0,                                          //tp_as_number
        0,                                          //tp_as_sequence
        0,                                          //tp_as_mapping
        0,                                          //tp_hash 
        0,                                          //tp_call
        0,                                          //tp_str
        0,                                          //tp_getattro
        0,                                          //tp_setattro
        0,                                          //tp_as_buffer
        Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   //tp_flags
        "Transform",                                //tp_doc 
        0,                                          //tp_traverse 
        0,                                          //tp_clear 
        0,                                          //tp_richcompare 
        0,                                          //tp_weaklistoffset 
        0,                                          //tp_iter 
        0,                                          //tp_iternext 
        PyOCIO_Transform_methods,                   //tp_methods 
        0,                                          //tp_members 
        0,                                          //tp_getset 
        0,                                          //tp_base 
        0,                                          //tp_dict 
        0,                                          //tp_descr_get 
        0,                                          //tp_descr_set 
        0,                                          //tp_dictoffset 
        (initproc) PyOCIO_Transform_init,           //tp_init 
        0,                                          //tp_alloc 
        0,                                          //tp_new 
        0,                                          //tp_free
        0,                                          //tp_is_gc
        0,                                          //tp_bases
        0,                                          //tp_mro
        0,                                          //tp_cache
        0,                                          //tp_subclasses
        0,                                          //tp_weaklist
        0,                                          //tp_del
        #if PY_VERSION_HEX > 0x02060000
        0,                                          //tp_version_tag
        #endif
    };
    
    ///////////////////////////////////////////////////////////////////////////
    ///
    
    namespace
    {
        ///////////////////////////////////////////////////////////////////////
        ///
        int PyOCIO_Transform_init( PyOCIO_Transform *self, PyObject * /*args*/, PyObject * /*kwds*/ )
        {
            ///////////////////////////////////////////////////////////////////
            /// init pyobject fields
            
            self->constcppobj = new ConstTransformRcPtr();
            self->cppobj = new TransformRcPtr();
            self->isconst = true;
            
            std::string message = "Base Transforms class can not be instantiated.";
            PyErr_SetString( PyExc_RuntimeError, message.c_str() );
            return -1;
        }
        
        ////////////////////////////////////////////////////////////////////////
        
        void PyOCIO_Transform_delete( PyOCIO_Transform *self, PyObject * /*args*/ )
        {
            delete self->constcppobj;
            delete self->cppobj;
            
            self->ob_type->tp_free((PyObject*)self);
        }
        
        ////////////////////////////////////////////////////////////////////////
        
        PyObject * PyOCIO_Transform_isEditable( PyObject * self )
        {
            return PyBool_FromLong(IsPyTransformEditable(self));
        }
        
        PyObject * PyOCIO_Transform_createEditableCopy( PyObject * self )
        {
            try
            {
                ConstTransformRcPtr transform = GetConstTransform(self, true);
                TransformRcPtr copy = transform->createEditableCopy();
                
                PyOCIO_Transform * pycopy = PyTransform_New(copy);
                pycopy->constcppobj = new ConstTransformRcPtr();
                pycopy->cppobj = new TransformRcPtr();
                *pycopy->cppobj = copy;
                pycopy->isconst = false;
                
                return (PyObject *) pycopy;
            }
            catch(...)
            {
                Python_Handle_Exception();
                return NULL;
            }
        }
        
        ////////////////////////////////////////////////////////////////////////
        
        
        PyObject * PyOCIO_Transform_getDirection( PyObject * self )
        {
            try
            {
                ConstTransformRcPtr transform = GetConstTransform(self, true);
                TransformDirection dir = transform->getDirection();
                return PyString_FromString( TransformDirectionToString( dir ) );
            }
            catch(...)
            {
                Python_Handle_Exception();
                return NULL;
            }
        }
        
        PyObject * PyOCIO_Transform_setDirection( PyObject * self, PyObject * args )
        {
            try
            {
                TransformDirection dir;
                if (!PyArg_ParseTuple(args,"O&:setDirection",
                    ConvertPyObjectToTransformDirection, &dir)) return NULL;
                
                TransformRcPtr transform = GetEditableTransform(self);
                transform->setDirection( dir );
                
                Py_RETURN_NONE;
            }
            catch(...)
            {
                Python_Handle_Exception();
                return NULL;
            }
        }
    }

}
OCIO_NAMESPACE_EXIT
