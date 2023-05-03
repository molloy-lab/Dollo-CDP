#include <boost/dynamic_bitset.hpp>



class Bipartition {

    public:
        Bipartition() {
        }

        Bipartition(boost::dynamic_bitset<> bitset){
            bip = bitset;
        }

        Bipartition(std::map<std::string, unsigned int> &label2index) {
            bip = boost::dynamic_bitset<>(label2index.size());
        }

        Bipartition(std::map<std::string, unsigned int> &label2index, std::vector<std::string> &subset) {
            bip = boost::dynamic_bitset<>(label2index.size());
            for (std::string label : subset) {
                bip.set(label2index[label]);
            }
        }

        Bipartition(std::vector<unsigned int> &subset, unsigned int size) {
            bip = boost::dynamic_bitset<>(size);
            for (unsigned int idx : subset) {
                bip.set(idx);
            }
        }

        Bipartition complement(const Bipartition &a) const {
            return Bipartition(a.bip - bip);
        }

        Bipartition other_child(const Bipartition &a) const {
            return Bipartition(a.bip ^ bip);
        }

        ~Bipartition() {
        }

        boost::dynamic_bitset<> get_bitset() const {
            return bip;
        }
        
        void set_bitset(boost::dynamic_bitset<> ib) {
            bip = ib;
        }

        std::string to_string() const {
            std::string s;
            boost::to_string(bip, s);
            return s;
        }
        std::string to_labels(std::vector<std::string> labels) const {
            std::string s;
            for (int i = 0; i < labels.size(); i ++) {
                if (bip.test(i)) s += " " + labels[i]; 
            }
            return s.substr(1, std::string::npos);
        }

        bool equivalent(const Bipartition &a, const Bipartition &b) {
            if (bip == b.bip) return true;
            if (a.bip == (bip ^ b.bip)) return true;
            return false;
        }
  /*
   index_t count() const {
            return bip.count();
        }
  */
        bool contain(std::map<std::string, unsigned int> &label2index, std::string &label) const {
            //std::cout << "bipartion class maps: "<< label << "-> " << label2index[label] << std::endl;
            return bip.test(label2index[label]);
        }

        int size() {
            return bip.size();
        }   

        int count() {
            return bip.count();
        }   

        int any() {
            return bip.any();
        }   
        void flip() {
            bip = ~bip;
        }


        inline bool contain_index(unsigned int ind) const {
            return bip.test(ind);
        }

    private:
        boost::dynamic_bitset<> bip;
};

bool operator==(const Bipartition &a, const Bipartition &b) {
    return a.get_bitset() == b.get_bitset();
}

namespace boost {
  std::size_t hash_value(Bipartition const& bip) {
    return std::hash<std::string>()(bip.to_string());
  }
}
template<>
struct std::hash<Bipartition> {
  unsigned int operator()(const Bipartition &bip) const {
    return std::hash<std::string>()(bip.to_string());
  }
};

