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

#include <LibParticle/ParticleSerialization.h>
#include <LibParticle/ParticleHelpers.h>

#include <map>
#include <algorithm>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
void ParticleSerialization::setFixedAttribute(const String& attrName, T value) {
    NT_REQUIRE(m_FixedAttributes.find(attrName) != m_FixedAttributes.end());
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    attr->buffer.setData(value);
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
void ParticleSerialization::setFixedAttribute(const String& attrName, const StdVT<T>& values) {
    NT_REQUIRE(m_FixedAttributes.find(attrName) != m_FixedAttributes.end());
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(values.size() == static_cast<size_t>(attr->count) && sizeof(T) == attr->typeSize());
    attr->buffer.setData(values, false);
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
void ParticleSerialization::setFixedAttribute(const String& attrName, T* values) {
    NT_REQUIRE(m_FixedAttributes.find(attrName) != m_FixedAttributes.end());
    auto& attr = m_FixedAttributes[attrName];
    attr->buffer.setData((const void*)values, attr->count * attr->typeSize());
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void ParticleSerialization::setFixedAttribute(const String& attrName, const VecX<N, T>& value) {
    NT_REQUIRE(m_FixedAttributes.find(attrName) != m_FixedAttributes.end());
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(static_cast<UInt>(N) == attr->count && sizeof(T) == attr->typeSize());
    attr->buffer.setData((const void*)glm::value_ptr(value), sizeof(T) * N);
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void ParticleSerialization::setFixedAttribute(const String& attrName, const MatXxX<N, T>& value) {
    NT_REQUIRE(m_FixedAttributes.find(attrName) != m_FixedAttributes.end());
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(static_cast<UInt>(N * N) == attr->count && sizeof(T) == attr->typeSize());
    attr->buffer.setData((const void*)glm::value_ptr(value), sizeof(T) * N * N);
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
void ParticleSerialization::setParticleAttribute(const String& attrName, const StdVT<T>& values) {
    NT_REQUIRE(m_ParticleAttributes.find(attrName) != m_ParticleAttributes.end() && values.size() == static_cast<size_t>(m_nParticles));
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(static_cast<UInt>(values.size()) == m_nParticles * attr->count);
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            NT_REQUIRE(sizeof(T) == attr->typeSize());
            attr->buffer.setData(values, false);
        } else {
            ParticleHelpers::compress(values, attr->buffer, false);
        }
    } else {
        NT_REQUIRE(sizeof(T) == attr->typeSize());
        attr->buffer.setData(values, false);
    }
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
void ParticleSerialization::setParticleAttribute(const String& attrName, const StdVT<StdVT<T>>& values) {
    NT_REQUIRE(m_ParticleAttributes.find(attrName) != m_ParticleAttributes.end() && values.size() == static_cast<size_t>(m_nParticles));
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(attr->type == TypeVectorInt || attr->type == TypeVectorUInt || attr->type == TypeVectorReal);

    NT_REQUIRE(sizeof(T) == attr->typeSize());
    attr->buffer.setData(values, false);
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void ParticleSerialization::setParticleAttribute(const String& attrName, const StdVT<VecX<N, T>>& values) {
    NT_REQUIRE(m_ParticleAttributes.find(attrName) != m_ParticleAttributes.end() && values.size() == static_cast<size_t>(m_nParticles));
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(static_cast<UInt>(N) == attr->count);
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            NT_REQUIRE(sizeof(T) == attr->typeSize());
            attr->buffer.setData(values, false);
        } else {
            ParticleHelpers::compress(values, attr->buffer, false);
        }
    } else {
        NT_REQUIRE(sizeof(T) == attr->typeSize());
        attr->buffer.setData(values, false);
    }
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void ParticleSerialization::setParticleAttribute(const String& attrName, const StdVT<MatXxX<N, T>>& values) {
    NT_REQUIRE(m_ParticleAttributes.find(attrName) != m_ParticleAttributes.end() && values.size() == static_cast<size_t>(m_nParticles));
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(static_cast<UInt>(N * N) == attr->count);
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            NT_REQUIRE(sizeof(T) == attr->typeSize());
            attr->buffer.setData((const void*)values.data(), values.size() * sizeof(T) * N * N);
        } else {
            ParticleHelpers::compress(values, attr->buffer, false);
        }
    } else {
        NT_REQUIRE(sizeof(T) == attr->typeSize());
        attr->buffer.setData((const void*)values.data(), values.size() * sizeof(T) * N * N);
    }
    attr->bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getFixedAttribute(const String& attrName, T& value) {
    if(m_FixedAttributes.find(attrName) == m_FixedAttributes.end()) {
        return false;
    }
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    attr->buffer.getData(value);
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getFixedAttribute(const String& attrName, T* values) {
    if(m_FixedAttributes.find(attrName) == m_FixedAttributes.end()) {
        return false;
    }
    auto& attr = m_FixedAttributes[attrName];
    attr->buffer.getData((void*)values, attr->count * attr->typeSize());
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getFixedAttribute(const String& attrName, StdVT<T>& values) {
    if(m_FixedAttributes.find(attrName) == m_FixedAttributes.end()) {
        return false;
    }
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    values.resize(attr->count);
    attr->buffer.getData(values, 0, attr->count);
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool ParticleSerialization::getFixedAttribute(const String& attrName, VecX<N, T>& value) {
    if(m_FixedAttributes.find(attrName) == m_FixedAttributes.end()) {
        return false;
    }
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(N == attr->count && sizeof(T) == attr->typeSize());
    attr->buffer.getData((void*)glm::value_ptr(value), sizeof(T) * N);
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool ParticleSerialization::getFixedAttribute(const String& attrName, MatXxX<N, T>& value) {
    if(m_FixedAttributes.find(attrName) == m_FixedAttributes.end()) {
        return false;
    }
    auto& attr = m_FixedAttributes[attrName];
    NT_REQUIRE(N * N == attr->count && sizeof(T) == attr->typeSize());
    attr->buffer.getData((void*)glm::value_ptr(value), sizeof(T) * N * N);
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getParticleAttribute(const String& attrName, T* values) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            attr->buffer.getData((void*)values, m_nParticles * attr->typeSize() * attr->count);
        } else {
            NT_REQUIRE(sizeof(T) == attr->typeSize());
            if(attr->count == 1) {
                StdVT<T> tmp;
                ParticleHelpers::decompress(tmp, attr->buffer, m_nParticles);
                memcpy(values, tmp.data(), attr->typeSize() * tmp.size());
            } else {
                if(attr->count == 2) {
                    StdVT<Vec2<T>> tmp;
                    ParticleHelpers::decompress(tmp, attr->buffer, m_nParticles);
                    memcpy(values, tmp.data(), attr->typeSize() * tmp.size());
                } else {
                    if(attr->count == 3) {
                        StdVT<Vec3<T>> tmp;
                        ParticleHelpers::decompress(tmp, attr->buffer, m_nParticles);
                        memcpy(values, tmp.data(), attr->typeSize() * tmp.size());
                    } else {
                        StdVT<Vec4<T>> tmp;
                        ParticleHelpers::decompress(tmp, attr->buffer, m_nParticles);
                        memcpy(values, tmp.data(), m_nParticles * attr->typeSize() * attr->count);
                    }
                }
            }
        }
    } else {
        attr->buffer.getData((void*)values, m_nParticles * attr->typeSize() * attr->count);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getParticleAttribute(const String& attrName, StdVT<T>& values) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    values.resize(m_nParticles);

    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            attr->buffer.getData(values, 0, m_nParticles);
        } else {
            ParticleHelpers::decompress(values, attr->buffer, m_nParticles);
        }
    } else {
        attr->buffer.getData(values, 0, m_nParticles);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getParticleAttribute(const String& attrName, StdVT<StdVT<T>>& values) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    values.resize(m_nParticles);
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            attr->buffer.getData(values, 0, m_nParticles);
        } else {
            ParticleHelpers::decompress(values, attr->buffer, m_nParticles);
        }
    } else {
        attr->buffer.getData(values, 0, m_nParticles);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool ParticleSerialization::getParticleAttribute(const String& attrName, StdVT<VecX<N, T>>& values) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(N == attr->count && sizeof(T) == attr->typeSize());
    values.resize(m_nParticles);
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            attr->buffer.getData(values, 0, m_nParticles);
        } else {
            ParticleHelpers::decompress(values, attr->buffer, m_nParticles);
        }
    } else {
        attr->buffer.getData(values, 0, m_nParticles);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool ParticleSerialization::getParticleAttribute(const String& attrName, StdVT<MatXxX<N, T>>& values) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(N * N == attr->count && sizeof(T) == attr->typeSize());
    values.resize(m_nParticles);
    if constexpr(std::is_floating_point_v<T>) {
        if(attr->type != TypeCompressedReal) {
            attr->buffer.getData((void*)values.data(), m_nParticles * sizeof(T) * N * N);
        } else {
            ParticleHelpers::decompress(values, attr->buffer, m_nParticles);
        }
    } else {
        attr->buffer.getData((void*)values.data(), m_nParticles * sizeof(T) * N * N);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getParticleAttributeCompressed(const String& attrName, StdVT_UInt16& values, T& dMin, T& dMax) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    NT_REQUIRE(attr->type == TypeCompressedReal);
    values.resize(m_nParticles);

    ////////////////////////////////////////////////////////////////////////////////
    float  dMinf, dMaxf;
    UInt64 segmentStart = 0;
    UInt64 segmentSize;

    segmentSize = sizeof(float);
    memcpy(&dMinf, &attr->buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;
    memcpy(&dMaxf, &attr->buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;

    segmentSize = m_nParticles * sizeof(UInt16) * attr->count;
    NT_REQUIRE(segmentStart + segmentSize == attr->buffer.size());
    memcpy(values.data(), &attr->buffer.data()[segmentStart], segmentSize);

    dMin = static_cast<T>(dMinf);
    dMax = static_cast<T>(dMaxf);
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool ParticleSerialization::getParticleAttributeCompressed(const String& attrName, StdVT_UInt16& values, VecX<N, T>& dMin, VecX<N, T>& dMax) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(N == attr->count && sizeof(T) == attr->typeSize());
    NT_REQUIRE(attr->type == TypeCompressedReal);
    values.resize(m_nParticles * N);

    ////////////////////////////////////////////////////////////////////////////////
    VecX<N, float> dMinf, dMaxf;
    UInt64         segmentStart = 0;
    UInt64         segmentSize  = sizeof(float) * N;

    memcpy(&dMinf, &attr->buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;
    memcpy(&dMaxf, &attr->buffer.data()[segmentStart], segmentSize);
    segmentStart += segmentSize;

    segmentSize = m_nParticles * sizeof(UInt16) * N;
    NT_REQUIRE(segmentStart + segmentSize == attr->buffer.size());
    memcpy(values.data(), &attr->buffer.data()[segmentStart], segmentSize);

    for(Int d = 0; d < N; ++d) {
        dMin[d] = static_cast<T>(dMinf[d]);
        dMax[d] = static_cast<T>(dMaxf[d]);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool ParticleSerialization::getParticleAttributeCompressed(const String& attrName, StdVT<StdVT_UInt16>& values, StdVT<T>& dMin, StdVT<T>& dMax) {
    if(m_ParticleAttributes.find(attrName) == m_ParticleAttributes.end()) {
        return false;
    }
    auto& attr = m_ParticleAttributes[attrName];
    NT_REQUIRE(sizeof(T) == attr->typeSize());
    NT_REQUIRE(attr->type == TypeCompressedReal);
    values.resize(m_nParticles);
    dMin.resize(m_nParticles);
    dMax.resize(m_nParticles);

    ////////////////////////////////////////////////////////////////////////////////
    float  dMinf;
    float  dMaxf;
    UInt64 segmentStart = 0;
    UInt64 segmentSize;

    for(UInt i = 0; i < m_nParticles; ++i) {
        UInt iSize;
        segmentSize = sizeof(UInt);
        memcpy(&iSize, &attr->buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;

        segmentSize = sizeof(float);
        memcpy(&dMinf, &attr->buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;
        memcpy(&dMaxf, &attr->buffer.data()[segmentStart], segmentSize);
        segmentStart += segmentSize;

        dMin[i] = static_cast<T>(dMinf);
        dMax[i] = static_cast<T>(dMaxf);

        values[i].resize(iSize);
        segmentSize = iSize * sizeof(UInt16);
        NT_REQUIRE(segmentStart + segmentSize <= attr->buffer.size());
        memcpy(values[i].data(), &attr->buffer.data()[segmentStart], segmentSize);
    }

    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void ParticleSerialization::saveParticle(const String& fileName, const StdVT<VecX<N, T>>& positions, T particleRadius, bool bCompress /*= true*/) {
    ParticleSerialization particleWriter;
    particleWriter.addFixedAttribute<T>("particle_radius", ParticleSerialization::TypeReal, 1);
    if(bCompress) {
        particleWriter.addParticleAttribute<T>("particle_position", ParticleSerialization::TypeCompressedReal, 3);
    } else {
        particleWriter.addParticleAttribute<T>("particle_position", ParticleSerialization::TypeReal, 3);
    }
    particleWriter.setNParticles(positions.size());
    particleWriter.setFixedAttribute("particle_radius", particleRadius);
    particleWriter.setParticleAttribute("particle_position", positions);
    particleWriter.flushAsync(fileName);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool ParticleSerialization::loadParticle(const String& fileName, StdVT<VecX<N, T>>& positions, T& particleRadius) {
    ParticleSerialization particleReader;
    if(!particleReader.read(fileName)) {
        return false;
    }

    T tmpRadius;
    NT_REQUIRE(particleReader.getFixedAttribute("particle_radius", tmpRadius));
    if(particleRadius > 0) {
        if(std::abs(tmpRadius - particleRadius) > MEpsilon<T>()) {
            return false;
        }
    } else {
        particleRadius = tmpRadius;
    }
    NT_REQUIRE(particleReader.getParticleAttribute("particle_position", positions));

    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// non-template functions
String ParticleSerialization::Attribute::typeName() {
    switch(type) {
        case TypeChar:
            return String("char");
        case TypeUInt16:
            return String("uint16");
        case TypeInt:
            return String("int");
        case TypeUInt:
            return String("uint");
        case TypeReal:
            return String("real");
        case TypeCompressedReal:
            return String("compressed_real");
        case TypeVectorChar:
            return String("vector_char");
        case TypeVectorInt:
            return String("vector_int");
        case TypeVectorUInt:
            return String("vector_uint");
        case TypeVectorReal:
            return String("vector_real");
        case TypeVectorCompressedReal:
            return String("vector_compressed_real");
        default:
            NT_DIE_UNKNOWN_ERROR
    }

    return String("");
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
size_t ParticleSerialization::Attribute::typeSize() {
    return static_cast<size_t>(size);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void ParticleSerialization::clearData() {
    m_nParticles = 0;

    for(auto& attr : m_FixedAttributes) {
        attr.second->buffer.clearBuffer();
        attr.second->bReady = false;
    }
    for(auto& attr : m_ParticleAttributes) {
        attr.second->buffer.clearBuffer();
        attr.second->bReady = false;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void ParticleSerialization::flushAsync(Int fileID) {
    NT_REQUIRE(m_DataIO != nullptr);
    m_DataIO->createOutputFolders();
    flushAsync(m_DataIO->getFilePath(fileID));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void ParticleSerialization::flushAsync(const String& fileName) {
    NT_REQUIRE(m_nParticles > 0 && m_FixedAttributes.size() > 0);
    if(m_Logger != nullptr) {
        buildAttrNameList();
        String str = String("Saving file: "); str += fileName;
        str += String(" ("); str += Formatters::toString(static_cast<double>(computeBufferSize()) / 1048576.0); str += String(" MBs)");
        m_Logger->printLog(str);
        str = String("File data: "); str += m_AttributeNameList;
        m_Logger->printLogIndent(str);
    }

    waitForBuffers();
    m_WriteFutureObj = std::async(std::launch::async, [&, fileName]() {
                                      std::ofstream opf(fileName, std::ios::binary | std::ios::out);
                                      if(!opf.is_open()) {
                                          if(m_Logger != nullptr) {
                                              m_Logger->printError("Cannot write file: " + fileName);
                                          }
                                          return;
                                      }

                                      ////////////////////////////////////////////////////////////////////////////////
                                      writeHeader(opf);
                                      for(auto& kv : m_FixedAttributes) {
                                          NT_REQUIRE(kv.second->bReady);
                                          opf.write((char*)kv.second->buffer.data(), kv.second->buffer.size());
                                      }
                                      for(auto& kv : m_ParticleAttributes) {
                                          NT_REQUIRE(kv.second->bReady);
                                          opf.write((char*)kv.second->buffer.data(), kv.second->buffer.size());
                                      }
                                      opf.close();
                                  });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
size_t ParticleSerialization::computeBufferSize() {
    size_t totalSize = 0;
    for(auto& kv : m_FixedAttributes) {
        totalSize += kv.second->buffer.size();
    }
    for(auto& kv : m_ParticleAttributes) {
        totalSize += kv.second->buffer.size();
    }

    return totalSize;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void ParticleSerialization::buildAttrNameList() {
    if(m_AttributeNameList.empty()) {
        for(auto& kv : m_FixedAttributes) {
            NT_REQUIRE_MSG(kv.second->bReady, kv.first + String(" attribute is not set!"));
            m_AttributeNameList += kv.first;
            m_AttributeNameList += String(", ");
        }
        for(auto& kv : m_ParticleAttributes) {
            NT_REQUIRE_MSG(kv.second->bReady, kv.first + String(" attribute is not set!"));
            m_AttributeNameList += kv.first;
            m_AttributeNameList += String(", ");
        }
        ////////////////////////////////////////////////////////////////////////////////
        m_AttributeNameList.erase(m_AttributeNameList.find_last_of(","), m_AttributeNameList.size()); // remove last ',' character
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void ParticleSerialization::writeHeader(std::ofstream& opf) {
    const std::locale& fixLoc = std::locale("C");
    opf.imbue(fixLoc);

    opf << "BananaParticleData\n";
    for(auto& kv : m_FixedAttributes) {
        if(kv.second->bReady) {
            opf << "FixedAttribute " << kv.second->name << " " << kv.second->typeName() << " " << kv.second->typeSize() << " " << kv.second->count << " " << kv.second->buffer.size() << "\n";
        }
    }

    for(auto& kv : m_ParticleAttributes) {
        if(kv.second->bReady) {
            opf << "ParticleAttribute " << kv.second->name << " " << kv.second->typeName() << " " << kv.second->typeSize() << " " << kv.second->count << " " << kv.second->buffer.size() << "\n";
        }
    }

    opf << "NParticles " << m_nParticles << "\n";
    opf << "EndHeader.\n";
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool ParticleSerialization::read(Int fileID, const StdVT<String>& readAttributes /*= {}*/, bool bStopIfFailed /*= true*/) {
    NT_REQUIRE(m_DataIO != nullptr);
    const String fileName = m_DataIO->getFilePath(fileID);
    return read(fileName, readAttributes, bStopIfFailed);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool ParticleSerialization::read(const String& fileName, const StdVT<String>& readAttributes /*= {}*/, bool bStopIfFailed /*= true*/) {
    std::ifstream ipf(fileName, std::ios::binary | std::ios::in);
    if(!ipf.is_open()) {
        if(m_Logger != nullptr) {
            m_Logger->printError("Cannot read file: " + fileName);
        }
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    clearData();
    m_ByteRead = 0;
    m_ReadAttributeDataSizeMap.clear();
    m_bReadAttributeMap.clear();
    if(!readHeader(ipf)) {
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    if(readAttributes.size() > 0) {
        // set all to false
        for(auto& kv : m_bReadAttributeMap) {
            kv.second = false;
        }

        // only given attributes will be true
        for(auto& attrName : readAttributes) {
            m_bReadAttributeMap[attrName] = true;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    size_t cursor = ipf.tellg();
    for(auto& kv : m_FixedAttributes) {
        if(m_bReadAttributeMap[kv.second->name]) {
            bool success = readAttribute(kv.second, ipf, cursor);
            if(!success && bStopIfFailed) {
                return false;
            }
            cursor = ipf.tellg();
        } else {
            size_t attrDataSize = m_ReadAttributeDataSizeMap[kv.second->name];
            cursor += attrDataSize;
        }
    }

    for(auto& kv : m_ParticleAttributes) {
        if(m_bReadAttributeMap[kv.second->name]) {
            bool success = readAttribute(kv.second, ipf, cursor);
            if(!success && bStopIfFailed) {
                return false;
            }
            cursor = ipf.tellg();
        } else {
            size_t attrDataSize = m_ReadAttributeDataSizeMap[kv.second->name];
            cursor += attrDataSize;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    ipf.close();
    if(m_Logger != nullptr) {
        String str = String("Read file: "); str += fileName;
        str += String(" ("); str += Formatters::toString(static_cast<double>(m_ByteRead) / 1048576.0); str += String(" MBs)");
        m_Logger->printLog(str);
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool ParticleSerialization::readHeader(Int fileID, const StdVT<String>& readAttributes /*= {}*/, bool bStopIfFailed /*= true*/) {
    NT_REQUIRE(m_DataIO != nullptr);
    const String fileName = m_DataIO->getFilePath(fileID);
    return readHeader(fileName, readAttributes, bStopIfFailed);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool ParticleSerialization::readHeader(const String& fileName, const StdVT<String>& readAttributes /*= {}*/, bool bStopIfFailed /*= true*/) {
    std::ifstream ipf(fileName, std::ios::binary | std::ios::in);
    if(!ipf.is_open()) {
        if(m_Logger != nullptr) {
            m_Logger->printError("Cannot read file: " + fileName);
        }
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    clearData();
    m_ByteRead = 0;
    m_ReadAttributeDataSizeMap.clear();
    m_bReadAttributeMap.clear();
    if(!readHeader(ipf)) {
        return false;
    }
    //todo: this function is not complte
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool ParticleSerialization::readHeader(std::ifstream& ipf) {
    String line;
    bool   gotMagic = false;

    auto getType = [&](const String& typeName) -> DataType {
                       if(typeName == "char") {
                           return TypeChar;
                       }
                       if(typeName == "uint16") {
                           return TypeUInt16;
                       }
                       if(typeName == "int") {
                           return TypeInt;
                       }
                       if(typeName == "uint") {
                           return TypeUInt;
                       }
                       if(typeName == "real") {
                           return TypeReal;
                       }
                       if(typeName == "compressed_real") {
                           return TypeCompressedReal;
                       }
                       if(typeName == "vector_char") {
                           return TypeVectorChar;
                       }
                       if(typeName == "vector_int") {
                           return TypeVectorInt;
                       }
                       if(typeName == "vector_uint") {
                           return TypeVectorUInt;
                       }
                       if(typeName == "vector_real") {
                           return TypeVectorReal;
                       }
                       if(typeName == "vector_compressed_real") {
                           return TypeVectorCompressedReal;
                       }
                       return TypeInt;
                   };

    while(std::getline(ipf, line)) {
        std::istringstream ls(line);
        String             token;
        ls >> token;

        if(token == "BananaParticleData") {
            gotMagic = true;
        } else if(token == "EndHeader.") {
            break;
        } else if(token == "NParticles") {
            ls >> m_nParticles;
        } else if(token == "FixedAttribute" || token == "ParticleAttribute") {
            String attrName, typeName;
            Int    typeSize;
            UInt   count;
            size_t dataSize;
            ls >> attrName >> typeName;
            ls >> typeSize >> count;
            ls >> dataSize;
            m_ReadAttributeDataSizeMap[attrName] = dataSize;
            m_bReadAttributeMap[attrName]        = true;

            if(token == "FixedAttribute") {
                m_FixedAttributes[attrName] = std::make_shared<Attribute>(attrName, getType(typeName), static_cast<ElementSize>(typeSize), count);
            } else {
                m_ParticleAttributes[attrName] = std::make_shared<Attribute>(attrName, getType(typeName), static_cast<ElementSize>(typeSize), count);
            }
        } else {
            return false;
        }
    }

    m_ByteRead = ipf.tellg();
    return gotMagic;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool ParticleSerialization::readAttribute(SharedPtr<Attribute>& attr, std::ifstream& ipf, size_t cursor) {
    if(m_ReadAttributeDataSizeMap.find(attr->name) == m_ReadAttributeDataSizeMap.end()) {
        return false;
    }
    size_t dataSize = m_ReadAttributeDataSizeMap[attr->name];
    ipf.seekg(cursor);
    attr->buffer.resize(dataSize);
    ipf.read((char*)attr->buffer.data(), dataSize);
    m_ByteRead += dataSize;

    return ipf.good();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __BNN_INSTANTIATE_SAVE_LOAD_PARTICLE(dimemsion, type)                                                                                \
    template void ParticleSerialization::saveParticle<dimemsion, type>(const String& fileName, const StdVT<VecX<dimemsion, type>>&positions, \
                                                                       type particleRadius, bool bCompress /*= true*/);                      \
    template bool ParticleSerialization::loadParticle<dimemsion, type>(const String& fileName, StdVT<VecX<dimemsion, type>>&positions,       \
                                                                       type & particleRadius);

__BNN_INSTANTIATE_SAVE_LOAD_PARTICLE(2, float)
__BNN_INSTANTIATE_SAVE_LOAD_PARTICLE(3, float)
__BNN_INSTANTIATE_SAVE_LOAD_PARTICLE(2, double)
__BNN_INSTANTIATE_SAVE_LOAD_PARTICLE(3, double)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(type)                                                                 \
    template void ParticleSerialization::setFixedAttribute<type>(const String& attrName, type value);               \
    template void ParticleSerialization::setFixedAttribute<type>(const String& attrName, type * value);             \
    template void ParticleSerialization::setFixedAttribute<const type>(const String& attrName, const type * value); \
    template void ParticleSerialization::setFixedAttribute<type>(const String& attrName, const StdVT<type>&value);

__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(float)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(double)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(Int)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(UInt)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(Int16)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(UInt16)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(Int64)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE(UInt64)

#define __BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(type)                                                        \
    template void ParticleSerialization::setFixedAttribute<2, type>(const String& attrName, const VecX<2, type>&value);   \
    template void ParticleSerialization::setFixedAttribute<3, type>(const String& attrName, const VecX<3, type>&value);   \
    template void ParticleSerialization::setFixedAttribute<4, type>(const String& attrName, const VecX<4, type>&value);   \
    template void ParticleSerialization::setFixedAttribute<2, type>(const String& attrName, const MatXxX<2, type>&value); \
    template void ParticleSerialization::setFixedAttribute<3, type>(const String& attrName, const MatXxX<3, type>&value); \
    template void ParticleSerialization::setFixedAttribute<4, type>(const String& attrName, const MatXxX<4, type>&value);

__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(double)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(Int)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(UInt)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(Int16)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(UInt16)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(Int64)
__BNN_INSTANTIATE_SET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(UInt64)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(type)                                                                \
    template void ParticleSerialization::setParticleAttribute<type>(const String& attrName, const StdVT<type>&value); \
    template void ParticleSerialization::setParticleAttribute<type>(const String& attrName, const StdVT<StdVT<type>>&value);

__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(float)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(double)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(Int)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(UInt)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(Int16)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(UInt16)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(Int64)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE(UInt64)

#define __BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(type)                                                               \
    template void ParticleSerialization::setParticleAttribute<2, type>(const String& attrName, const StdVT<VecX<2, type>>&value);   \
    template void ParticleSerialization::setParticleAttribute<3, type>(const String& attrName, const StdVT<VecX<3, type>>&value);   \
    template void ParticleSerialization::setParticleAttribute<4, type>(const String& attrName, const StdVT<VecX<4, type>>&value);   \
    template void ParticleSerialization::setParticleAttribute<2, type>(const String& attrName, const StdVT<MatXxX<2, type>>&value); \
    template void ParticleSerialization::setParticleAttribute<3, type>(const String& attrName, const StdVT<MatXxX<3, type>>&value); \
    template void ParticleSerialization::setParticleAttribute<4, type>(const String& attrName, const StdVT<MatXxX<4, type>>&value);

__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(double)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(Int)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(UInt)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(Int16)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(UInt16)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(Int64)
__BNN_INSTANTIATE_SET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(UInt64)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(type)                                                     \
    template bool ParticleSerialization::getFixedAttribute<type>(const String& attrName, type & value); \
    template bool ParticleSerialization::getFixedAttribute<type>(const String& attrName, type * value); \
    template bool ParticleSerialization::getFixedAttribute<type>(const String& attrName, StdVT<type>&value);

__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(float)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(double)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(Int)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(UInt)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(Int16)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(UInt16)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(Int64)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE(UInt64)

#define __BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(type)                                                  \
    template bool ParticleSerialization::getFixedAttribute<2, type>(const String& attrName, VecX<2, type>&value);   \
    template bool ParticleSerialization::getFixedAttribute<3, type>(const String& attrName, VecX<3, type>&value);   \
    template bool ParticleSerialization::getFixedAttribute<4, type>(const String& attrName, VecX<4, type>&value);   \
    template bool ParticleSerialization::getFixedAttribute<2, type>(const String& attrName, MatXxX<2, type>&value); \
    template bool ParticleSerialization::getFixedAttribute<3, type>(const String& attrName, MatXxX<3, type>&value); \
    template bool ParticleSerialization::getFixedAttribute<4, type>(const String& attrName, MatXxX<4, type>&value);

__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(double)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(Int)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(UInt)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(Int16)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(UInt16)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(Int64)
__BNN_INSTANTIATE_GET_FIXED_ATTRIBUTE_COMMON_VEC_DIM(UInt64)
////////////////////////////////////////////////////////////////////////////////
#define __BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(type)                                                          \
    template bool ParticleSerialization::getParticleAttribute<type>(const String& attrName, type * value);      \
    template bool ParticleSerialization::getParticleAttribute<type>(const String& attrName, StdVT<type>&value); \
    template bool ParticleSerialization::getParticleAttribute<type>(const String& attrName, StdVT<StdVT<type>>&value);

__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(float)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(double)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(Int)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(UInt)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(Int16)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(UInt16)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(Int64)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE(UInt64)

#define __BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(type)                                                         \
    template bool ParticleSerialization::getParticleAttribute<2, type>(const String& attr2ame, StdVT<VecX<2, type>>&value);   \
    template bool ParticleSerialization::getParticleAttribute<3, type>(const String& attr3ame, StdVT<VecX<3, type>>&value);   \
    template bool ParticleSerialization::getParticleAttribute<4, type>(const String& attr4ame, StdVT<VecX<4, type>>&value);   \
    template bool ParticleSerialization::getParticleAttribute<2, type>(const String& attr2ame, StdVT<MatXxX<2, type>>&value); \
    template bool ParticleSerialization::getParticleAttribute<3, type>(const String& attr3ame, StdVT<MatXxX<3, type>>&value); \
    template bool ParticleSerialization::getParticleAttribute<4, type>(const String& attr4ame, StdVT<MatXxX<4, type>>&value);

__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(float)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(double)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(Int)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(UInt)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(Int16)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(UInt16)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(Int64)
__BNN_INSTANTIATE_GET_PARTICLE_ATTRIBUTE_COMMON_VEC_DIM(UInt64)

//template bool ParticleSerialization::getParticleAttributeCompressed(const String& attrName, StdVT_UInt16& values, T& dMin, T& dMax);
//template bool ParticleSerialization::getParticleAttributeCompressed(const String& attrName, StdVT_UInt16& values, VecX<N, T>& dMin, VecX<N, T>& dMax);
//template bool ParticleSerialization::getParticleAttributeCompressed(const String& attrName, StdVT<StdVT_UInt16>& values, StdVT<T>& dMin, StdVT<T>& dMax);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
