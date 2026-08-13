//
// Created by Patrick Flege on 18/03/2026.
//

#ifndef KMERTYPE_H
#define KMERTYPE_H



enum class KmerType {
   CORE = 1, ACCESSORY = 2, UNIQUE = 3
};

inline std::ostream& operator<<(std::ostream& os, const KmerType& kmer_type) {
    switch (kmer_type) {
        case KmerType::CORE:
            os << "CORE" << std::endl;
            break;
        case KmerType::ACCESSORY:
            os << "ACCESSORY" << std::endl;
            break;
        case KmerType::UNIQUE:
            os << "UNIQUE" << std::endl;
            break;
    }
    return os;
};



#endif //KMERTYPE_H
