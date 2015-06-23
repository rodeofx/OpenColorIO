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

#include <tinyxml.h>
#include <stdio.h>
#include <vector>

#include <OpenColorIO/OpenColorIO.h>

#include "FileTransform.h"
#include "OpBuilders.h"
#include "MatrixOps.h"
#include "OpBuilders.h"
#include "CDLTransform.h"

OCIO_NAMESPACE_ENTER
{
    ////////////////////////////////////////////////////////////////
    
    namespace
    {
        /// An internal Op cache that has been read from a CTF file.
        class CachedOp
        {
        public:
            virtual void buildFinalOp(OpRcPtrVec &ops,
                                      const Config& config,
                                      TransformDirection dir) = 0;
        };

        typedef OCIO_SHARED_PTR<CachedOp> CachedOpRcPtr;
        typedef std::vector<CachedOpRcPtr> CachedOpRcPtrVec;

        ////////////////////////////////////////////////////////////////

        class XMLTagHandler;
        typedef OCIO_SHARED_PTR<XMLTagHandler> XMLTagHandlerRcPtr;

        class XMLTagHandler
        {
        public:
            static XMLTagHandlerRcPtr CreateHandlerForTagName(std::string text);

            virtual CachedOpRcPtr handleXMLTag(TiXmlElement * element) = 0;
        };

        /// To add support for a new CTF XML tag:
        ///
        /// 1. Copy the template FooTagHandler below and replace Foo with an
        /// appropriate name.
        ///
        /// 2. Modify the CreateHandlerForTagName factory method to return
        /// your new handler subclass when passed an appropriate XML tag string.

        /*
        class FooTagHandler : public XMLTagHandler
        {
            class FooCachedOp : public CachedOp
            {
            public:
                /// Data structures representing the operator go here.
                /// For example, the MatrixCachedOp contains just a float m_m44[16],
                /// to store the elements of the matrix.

                virtual void buildFinalOp(OpRcPtrVec &ops,
                                          const Config& config,
                                          TransformDirection dir) {
                    /// Turn the cached data into actual ops by passing the 
                    /// ops variable to a Create___Op function.

                    /// For example, the matrix operator just called
                    /// CreateMatrixOp(ops, &m_m44[0], dir);
                }
            };

            virtual CachedOpRcPtr handleXMLTag(TiXmlElement * element) {
                /// FooCachedOp * cachedOp = new FooCachedOp;

                /// Read the ASCII data from the Tiny XML element representation
                /// into the intermediate cached format.
                /// For example, the MatrixTagHandler just converts the ASCII
                /// values to floats and stores them in a MatrixCachedOp.               

                /// return CachedOpRcPtr(cachedOp);
            }
        };
        */

        /// Handles <matrix> tags in CTF files.
        class MatrixTagHandler : public XMLTagHandler
        {
            class MatrixCachedOp : public CachedOp
            {
            public:
                const unsigned short DIM_SIZE;
                const unsigned short MTX_SIZE;
                const unsigned short MTX_DIM;
                // Dimension and matrix members
                std::vector<unsigned short> m_dim;
                std::vector<float> m_m44;

                MatrixCachedOp() : DIM_SIZE(3), MTX_DIM(4), MTX_SIZE(16)
                {
                    // Default values
                    const unsigned short default_dimension[] = {3, 3, 3};
                    const float default_matrix[] = {1.f, 0.f, 0.f, 0.f,
                                                    0.f, 1.f, 0.f, 0.f,
                                                    0.f, 0.f, 1.f, 0.f,
                                                    0.f, 0.f, 0.f, 1.f};
                    // Init members to default values
                    this->m_dim.assign(&default_dimension[0],
                                       &default_dimension[0] + DIM_SIZE);
                    this->m_m44.assign(&default_matrix[0],
                                       &default_matrix[0] + MTX_SIZE);
                }

                /// Turns the cached data into an Op and adds it to the vector of
                /// Ops
                virtual void buildFinalOp(OpRcPtrVec &ops,
                                          const Config& config,
                                          TransformDirection dir) {
                    CreateMatrixOp(ops, &m_m44[0], dir);
                }

                // Parse string to fill vector considering input data size
                // Ex: fill the 4x4 matrix member with a 3x4 input matrix
                template <typename T>
                void fillVectorFromString(std::vector<T> &v, const char * str,
                                          const unsigned short inColSize)
                {
                    std::istringstream is(str);
                    int i = 0;
                    int j = 1;
                    int delta = this->MTX_DIM - inColSize;
                    while (!is.eof() && i < v.size()) {
                        T token;
                        is >> token;
                        v[i++] = token;
                        if (i >= inColSize*j + (j-delta))
                        {
                            i += delta;
                            ++j;
                        }
                    }
                }

                template <typename T>
                void fillVectorFromString(std::vector<T> &v, const char * str)
                {
                    fillVectorFromString(v, str, v.size());
                }

                // Check dimensions.
                // Specs possible dimensions 3x3 3, 3x4 3 ; OCIO matrix op can
                // handle 4x4 4.
                // TODO check if we really want to support 4x4 while not
                // defined by specs.
                void checkDimension()
                {
                    unsigned short valid_dims[3][3] = {{3, 3 ,3},
                                                       {3, 4, 3},
                                                       {4, 4, 4}};
                    for (int i=0; i < 3 ; ++i)
                    {
                        std::vector<unsigned short> valid_dim;
                        valid_dim.assign(&valid_dims[i][0],
                                         &valid_dims[i][0] + DIM_SIZE);
                        if (std::equal(this->m_dim.begin(),
                                       this->m_dim.begin() + DIM_SIZE,
                                       valid_dim.begin()))
                            return;
                    }
                    std::ostringstream os;
                    os << "Matrix Parsing Error: wrong dimension ("
                       << this->m_dim[0] << "x" << this->m_dim[1] << " "
                       << this->m_dim[2] <<  ")";
                    throw Exception(os.str().c_str());
                }
            };

            /// Parses the raw ASCII data and returns a CachedOp
            /// containing the data.
            virtual CachedOpRcPtr handleXMLTag(TiXmlElement * element) {
                MatrixCachedOp * cachedOp = new MatrixCachedOp;

                // Find the Array XML tag
                TiXmlElement *arrayElement = TiXmlHandle(element).
                                             FirstChildElement("Array").
                                             ToElement();

                if (!arrayElement){
                    std::ostringstream os;
                    os << "Matrix Parsing Error: could not find XML Array element !";
                    throw Exception(os.str().c_str());
                }

                // Read the matrix tokens into our matrix array
                // Get dim attribute
                cachedOp->fillVectorFromString(cachedOp->m_dim,
                                               arrayElement->Attribute("dim"),
                                               cachedOp->DIM_SIZE);
                cachedOp->checkDimension();
                // Get matrix values
                cachedOp->fillVectorFromString(cachedOp->m_m44,
                                               arrayElement->GetText(),
                                               cachedOp->m_dim[1]);

                return CachedOpRcPtr(cachedOp);
            }
        };

        /// Handles <ASC_CDL> tags
        // TODO Handle styles : Fwd, Rev, FwdNoClamp, RevNoClamp
        class CDLTagHandler : public XMLTagHandler
        {
            class CDLCachedOp : public CachedOp
            {
            public:
                CDLTransformRcPtr transform;

                CDLCachedOp ()
                {
                    transform = CDLTransform::Create();
                };

                virtual void buildFinalOp(OpRcPtrVec &ops,
                                          const Config& config,
                                          TransformDirection dir) {
                    BuildCDLOps(ops, config, *(this->transform), dir);
                }
            };

            virtual CachedOpRcPtr handleXMLTag(TiXmlElement * element) {
                CDLCachedOp * cachedOp = new CDLCachedOp;
                // Rename ASC_CDL into ColorCorrection
                element->SetValue("ColorCorrection");
                // Load XML data into the transform
                LoadCDL(cachedOp->transform.get(), element);
                return CachedOpRcPtr(cachedOp);
            }
        };

        /// A factory method to instantiate an appropriate XMLTagHandler for a given
        /// tag name. For example, passing "matrix" to this function should instantiate
        /// and return a MatrixTagHandler.
        XMLTagHandlerRcPtr XMLTagHandler::CreateHandlerForTagName(std::string text)
        {
            if (text.compare("Matrix") == 0) {
                return XMLTagHandlerRcPtr(new MatrixTagHandler());
            }else if (text.compare("ASC_CDL") == 0) {
                return XMLTagHandlerRcPtr(new CDLTagHandler());
            }
            return XMLTagHandlerRcPtr(static_cast<XMLTagHandler*>(NULL));
        }

        ////////////////////////////////////////////////////////////////

        class LocalCachedFile : public CachedFile
        {
        public:
            LocalCachedFile ()
            {

            };
            
            ~LocalCachedFile() {};
            
            CachedOpRcPtrVec m_cachedOps;

        };
        
        typedef OCIO_SHARED_PTR<LocalCachedFile> LocalCachedFileRcPtr;
        typedef OCIO_SHARED_PTR<TiXmlDocument> TiXmlDocumentRcPtr;
        
        
        class LocalFileFormat : public FileFormat
        {
        public:
            
            ~LocalFileFormat() {};
            
            virtual void GetFormatInfo(FormatInfoVec & formatInfoVec) const;
            
            virtual CachedFileRcPtr Read(std::istream & istream) const;
            
            virtual void BuildFileOps(OpRcPtrVec & ops,
                                      const Config& config,
                                      const ConstContextRcPtr & context,
                                      CachedFileRcPtr untypedCachedFile,
                                      const FileTransform& fileTransform,
                                      TransformDirection dir) const;
        };
        

        void LocalFileFormat::GetFormatInfo(FormatInfoVec & formatInfoVec) const
        {
            FormatInfo info;
            info.name = "Color Transform File";
            info.extension = "ctf";
            info.capabilities = FORMAT_CAPABILITY_READ;
            formatInfoVec.push_back(info);
        }
        
        // Try and load the format
        // Raise an exception if it can't be loaded.
        
        CachedFileRcPtr LocalFileFormat::Read(std::istream & istream) const
        {
            std::ostringstream rawdata;
            rawdata << istream.rdbuf();

            // We were asked to Read, so we should create a new cached file and fill it with
            // data from the .ctf file.
            LocalCachedFileRcPtr cachedFile = LocalCachedFileRcPtr(new LocalCachedFile());
            
            // Create a TinyXML representation of the raw data we were given
            TiXmlDocumentRcPtr doc = TiXmlDocumentRcPtr(new TiXmlDocument());
            doc->Parse(rawdata.str().c_str());
            
            if(doc->Error())
            {
                std::ostringstream os;
                os << "XML Parse Error. ";
                os << doc->ErrorDesc() << " (line ";
                os << doc->ErrorRow() << ", character ";
                os << doc->ErrorCol() << ")";
                throw Exception(os.str().c_str());
            }
            
            // Get the root element of the CTF file. It should be "ProcessList"
            // TODO: should probably verify this!
            TiXmlElement* rootElement = doc->RootElement();

            // Iterate through the children nodes
            TiXmlElement* currentElement = rootElement->FirstChildElement();
            while (currentElement) {
                std::string tagName = currentElement->Value();

                // Create an XMLTagHandler to handle this specific tag
                XMLTagHandlerRcPtr tagHandler = XMLTagHandler::CreateHandlerForTagName(tagName);
                if (tagHandler) {
                    CachedOpRcPtr cachedOp = tagHandler->handleXMLTag(currentElement);

                    // Store the CachedOp in our cached file, as we'll need it later
                    cachedFile->m_cachedOps.push_back(cachedOp);
                }
                
                // On to the next one
                currentElement = currentElement->NextSiblingElement();
            }
            
            return cachedFile;
        }
        
        void
        LocalFileFormat::BuildFileOps(OpRcPtrVec & ops,
                                      const Config& config,
                                      const ConstContextRcPtr & /*context*/,
                                      CachedFileRcPtr untypedCachedFile,
                                      const FileTransform& fileTransform,
                                      TransformDirection dir) const
        {
            LocalCachedFileRcPtr cachedFile = DynamicPtrCast<LocalCachedFile>(untypedCachedFile);
            
            // This should never happen.
            if(!cachedFile)
            {
                std::ostringstream os;
                os << "Cannot build .ctf Op. Invalid cache type.";
                throw Exception(os.str().c_str());
            }
            
            TransformDirection newDir = CombineTransformDirections(dir,
                fileTransform.getDirection());
            if(newDir == TRANSFORM_DIR_UNKNOWN)
            {
                std::ostringstream os;
                os << "Cannot build file format transform,";
                os << " unspecified transform direction.";
                throw Exception(os.str().c_str());
            }

            // Iterate through the cached file ops
            for (int i = 0; i < cachedFile->m_cachedOps.size(); i++) {
                CachedOpRcPtr thisOp = cachedFile->m_cachedOps[i];
                thisOp->buildFinalOp(ops, config, dir);
            }
        }
    }
    
    FileFormat * CreateFileFormatCTF()
    {
        return new LocalFileFormat();
    }
}
OCIO_NAMESPACE_EXIT

