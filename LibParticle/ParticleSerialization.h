//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTCodeBase                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Data/DataIO.h>
#include <LibCommon/Logger/Logger.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Banana Particle Data format =====>
////////////////////////////////////////////////////////////////////////////////
//
// "BananaParticleData"
// "<FixedAtributeName1>:" <DataType> <ElementSize> <Count> <DataSize>
//      ...
//
// "<ParticleAttribute1>:" <DataType> <ElementSize> <Count> <DataSize>
//      ...
//
// "NParticles:" <Num. of particles>
// "EndHeader."
// <DataOfParticleAttribute1>
// <DataOfParticleAttribute2>
// ...
//
////////////////////////////////////////////////////////////////////////////////
// =====> Banana Particle Data format
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class ParticleSerialization {
public:
    enum DataType {
        TypeChar,
        TypeUInt16,
        TypeInt,
        TypeUInt,
        TypeReal,
        TypeCompressedReal,
        TypeVectorChar,
        TypeVectorInt,
        TypeVectorUInt,
        TypeVectorReal,
        TypeVectorCompressedReal
    };

    enum ElementSize {
        Size8b  = sizeof(char),
        Size16b = sizeof(UInt16),
        Size32b = sizeof(float),
        Size64b = sizeof(double)
    };

    struct Attribute {
        DataBuffer  buffer;
        String      name;
        DataType    type;
        ElementSize size;
        UInt        count;
        bool        bReady;

        Attribute(const String& name_, DataType type_, ElementSize size_, UInt count_ = 1) :
            name(name_), type(type_), size(size_), count(count_), bReady(false) {}
        String typeName();
        size_t typeSize();
    };

public:
    ParticleSerialization() = default;
    ParticleSerialization(const String&            dataRootFolder,
                          const String&            dataFolder,
                          const String&            fileName,
                          const SharedPtr<Logger>& logger = nullptr) :
        m_DataIO(std::make_shared<DataIO>(dataRootFolder, dataFolder, fileName, String("bnn"), String("BananaParticleData"))), m_Logger(logger)
    {}

    virtual ~ParticleSerialization() { waitForBuffers(); }

    void setDataPath(const String& dataRootFolder,
                     const String& dataFolder,
                     const String& fileName) {
        m_DataIO = std::make_shared<DataIO>(dataRootFolder, dataFolder, fileName, String("bnn"), String("BananaParticleData"));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // functions for writing data
    template<class T>
    void addFixedAttribute(const String& attrName, DataType type, UInt count = 1) {
        addFixedAttribute(attrName, type, static_cast<ParticleSerialization::ElementSize>(sizeof(T)), count);
    }

    template<class T>
    void addParticleAttribute(const String& attrName, DataType type, UInt count = 1) {
        addParticleAttribute(attrName, type, static_cast<ParticleSerialization::ElementSize>(sizeof(T)), count);
    }

    void addFixedAttribute(const String& attrName, DataType type, ElementSize size, UInt count = 1) {
        NT_REQUIRE(type == TypeChar || type == TypeInt || type == TypeUInt || type == TypeReal);
        m_FixedAttributes[attrName] = std::make_shared<Attribute>(attrName, type, size, count);
    }

    void addParticleAttribute(const String& attrName, DataType type, ElementSize size, UInt count = 1) {
        m_ParticleAttributes[attrName] = std::make_shared<Attribute>(attrName, type, size, count);
    }

    ////////////////////////////////////////////////////////////////////////////////
    template<class IndexType> void setNParticles(IndexType nParticles) { m_nParticles = static_cast<UInt>(nParticles); }

    template<class T> void        setFixedAttribute(const String& attrName, T value);
    template<class T> void        setFixedAttribute(const String& attrName, T* values);
    template<class T> void        setFixedAttribute(const String& attrName, const StdVT<T>& values);
    template<Int N, class T> void setFixedAttribute(const String& attrName, const VecX<N, T>& value);
    template<Int N, class T> void setFixedAttribute(const String& attrName, const MatXxX<N, T>& value);

    template<class T> void        setParticleAttribute(const String& attrName, const StdVT<T>& values);
    template<class T> void        setParticleAttribute(const String& attrName, const StdVT<StdVT<T>>& values);
    template<Int N, class T> void setParticleAttribute(const String& attrName, const StdVT<VecX<N, T>>& values);
    template<Int N, class T> void setParticleAttribute(const String& attrName, const StdVT<MatXxX<N, T>>& values);

    void clearData();
    void flushAsync(Int fileID);
    void flushAsync(const String& fileName);
    void waitForBuffers() { if(m_WriteFutureObj.valid()) { m_WriteFutureObj.wait(); } }

    ////////////////////////////////////////////////////////////////////////////////
    // functions for reading data
    auto& getFixedAttributes() { return m_FixedAttributes; }
    const auto& getFixedAttributes() const { return m_FixedAttributes; }
    auto& getParticleAttributes() { return m_ParticleAttributes; }
    const auto& getParticleAttributes() const { return m_ParticleAttributes; }

    Int    getLatestFileIndex(Int maxIndex) const { return m_DataIO->getLatestFileIndex(maxIndex); }
    String getFilePath(Int fileID) { return m_DataIO->getFilePath(fileID); }
    bool read(Int fileID, const StdVT<String>& readAttributes = {}, bool bStopIfFailed = true);
    bool read(const String& fileName, const StdVT<String>& readAttributes       = {}, bool bStopIfFailed = true);
    bool readHeader(Int fileID, const StdVT<String>& readAttributes             = {}, bool bStopIfFailed = true);
    bool readHeader(const String& fileName, const StdVT<String>& readAttributes = {}, bool bStopIfFailed = true);
    size_t getBytesRead() const { return m_ByteRead; }
    UInt   nParticles() const { return m_nParticles; }

    bool hasFixedAttribute(const String& attrName) const { return m_FixedAttributes.find(attrName) != m_FixedAttributes.end(); }
    bool hasParticleAttribute(const String& attrName) const { return m_ParticleAttributes.find(attrName) != m_ParticleAttributes.end(); }
    bool hasAttribute(const String& attrName) const { return hasFixedAttribute(attrName) || hasParticleAttribute(attrName); }

    template<class T> bool        getFixedAttribute(const String& attrName, T& value);
    template<class T> bool        getFixedAttribute(const String& attrName, T* value);
    template<class T> bool        getFixedAttribute(const String& attrName, StdVT<T>& values);
    template<Int N, class T> bool getFixedAttribute(const String& attrName, VecX<N, T>& value);
    template<Int N, class T> bool getFixedAttribute(const String& attrName, MatXxX<N, T>& value);

    template<class T> bool        getParticleAttribute(const String& attrName, T* values);
    template<class T> bool        getParticleAttribute(const String& attrName, StdVT<T>& values);
    template<class T> bool        getParticleAttribute(const String& attrName, StdVT<StdVT<T>>& values);
    template<Int N, class T> bool getParticleAttribute(const String& attrName, StdVT<VecX<N, T>>& values);
    template<Int N, class T> bool getParticleAttribute(const String& attrName, StdVT<MatXxX<N, T>>& values);

    template<class T> bool        getParticleAttributeCompressed(const String& attrName, StdVT_UInt16& values, T& dMin, T& dMax);
    template<Int N, class T> bool getParticleAttributeCompressed(const String& attrName, StdVT_UInt16& values, VecX<N, T>& dMin, VecX<N, T>& dMax);
    template<class T> bool        getParticleAttributeCompressed(const String& attrName, StdVT<StdVT_UInt16>& values, StdVT<T>& dMin, StdVT<T>& dMax);

private:
    size_t computeBufferSize();
    void   buildAttrNameList();
    void   writeHeader(std::ofstream& opf);
    bool   readHeader(std::ifstream& ipf);
    bool   readAttribute(SharedPtr<Attribute>& attr, std::ifstream& ipf, size_t cursor);

    UInt m_nParticles;
    Map<String, SharedPtr<Attribute>> m_FixedAttributes;
    Map<String, SharedPtr<Attribute>> m_ParticleAttributes;
    String m_AttributeNameList;

    Map<String, size_t> m_ReadAttributeDataSizeMap;
    Map<String, bool>   m_bReadAttributeMap;
    size_t              m_ByteRead;

    SharedPtr<Logger> m_Logger;
    SharedPtr<DataIO> m_DataIO;
    std::future<void> m_WriteFutureObj;

public:
    template<Int N, class T> static void saveParticle(const String& fileName, const StdVT<VecX<N, T>>& positions, T particleRadius, bool bCompress = true);
    template<Int N, class T> static bool loadParticle(const String& fileName, StdVT<VecX<N, T>>& positions, T& particleRadius);
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
