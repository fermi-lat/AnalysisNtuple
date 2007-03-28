#include "ntupleWriterSvc/INTupleWriterSvc.h"
INTupleWriterSvc* rootTupleSvc;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class Item {
    public:
        Item(std::string name, char typecode=' ')
        {
            std::string type = rootTupleSvc->getItem(treename, name, m_pvalue);
            if( typecode==' ') {
                m_isFloat = type==rootType('F');
                if( !m_isFloat && type!=rootType('D')){
                    throw std::invalid_argument("McCoordsAlg: type of "+name+ " is not "+ rootType('F')+" or "+rootType('D'));
                }
            }else if( type!= rootType(typecode) ){
                throw std::invalid_argument("McCoordsAlg: type of "+name+ " is not "+ rootType(typecode));
            }
        }
        // Item behaves like a double
        operator double()const
        {
            return m_isFloat? *(float*)m_pvalue : *(double*)m_pvalue;
        }

        static std::string rootType(char code){
            if( code=='i') return "UInt_t";
            if( code=='I') return "Int_t";
            if( code=='F') return "Float_t";
            if( code=='D') return "Double_t";
            // todo: add more?
            return "unknown";
        }
        void* m_pvalue;
        bool m_isFloat;
    };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    template<typename T, char  typecode>
    class TypedItem : public Item {
    public:
        TypedItem(std::string name): Item(name, typecode){}
        T value() const{ return *static_cast<T*>(m_pvalue); }
        operator T()const{return value();}
    };
    template <typename T>
        void addItem(std::string name, const T & value)
    {
        rootTupleSvc->addItem(treename, name, &value);
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
