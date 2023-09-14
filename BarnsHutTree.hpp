#pragma once
#include <functional>
#include <iostream>
#include <cmath>
#include <limits>

namespace 
{
    template <class type, class... Args>
    type* Allocate(Args... args)
    {
        type* p = new (std::nothrow) type(args...);
        if(p==nullptr)
        {
            std::cout << "dont't have enough merty";
            std::exit(-1);
        }
        return p;
    }

    template <class type>
    type* Allocate_Array(unsigned int n)
    {
        type* p = new (std::nothrow) type[n];
        if(p==nullptr)
        {
            std::cout << "dont't have enough memty";
            std::exit(-1);
        }

        return p;
    }

    
    int ncell = 0;

    enum class Type : unsigned char
    {
        Body = 0,
        Cell,
    };

    struct Node
    {
        Type type;
        bool update;
        double* position;
        Node* next;
    };

    struct Body:public Node
    {
        unsigned int id;
    };
}


template <class posType, class AdditionalCellInfo>
class BarnesHutTree
{
public:
    struct Cell:public Node
    {
        bool flag = false;
        Node* more = nullptr;
        Node** subP;
        AdditionalCellInfo info;
        
        Cell()
        {

        }   

        Cell(unsigned char DIM):flag(true)
        {
            subP = new Node*[1<<DIM];
            position = new double[DIM];
        }

        void Construct(unsigned int DIM)
        {
            subP = new Node*[1<<DIM];
            position = new double[DIM];
            flag = true;
        }

        ~Cell()
        {
            if(flag)
            {   
                delete[] subP;
                delete[] position;
            }
        }

        void ClearPosition(unsigned char DIM)
        {
            for(int i = 0;i<DIM;++i)
                position[i] = 0;
        }
    };

private:
    Cell* root;//root pointer
    Body* bodies;
    Body* bodies2;
    Node* freeCell = nullptr;
    int NSUB;
    int actLen;
    Node** active;
    Cell* interact;
    double length;
    unsigned char DIM;
    std::function<void(int,int)> action;
    std::function<void(int,Cell*)> bc_action;
    std::function<bool(Cell* c, double psize, double* pmid)> accept;
    bool firstCall = true;

    void NewTree()
    {
        Node* p;
        if(!firstCall)
        {
            p = root;
            while(p != nullptr)
            {
                if(p->type == Type::Cell)
                {
                    p->next = freeCell;
                    freeCell = p;
                    p = static_cast<Cell*>(p)->more;
                }else
                {
                    p = p->next;
                }
            }
        }else
        {
            firstCall = false;
        }

        root = nullptr;
        ncell = 0;
    }

    Cell* MakeCell()
    {
        Cell* c;
        int i;

        if(freeCell == nullptr)
        {
            c = new Cell(DIM);
        }
        else
        {
            c = static_cast<Cell*>(freeCell);
            freeCell = c->next;
        }

        c->type = Type::Cell;
        c->update = false;
        for(int i = 0;i<NSUB;++i)
            c->subP[i] = nullptr;
        ++ncell;
        return c;
    }

    void ExpandBox()
    {
        rsize = 1;
        double dmax = 0,d;
        Body* p;
        for(p = bodies;p<bodies + num;++p)
        {
            for(int dim = 0;dim<DIM;++dim)
            {
                d = std::abs(p->position[dim]-root->position[dim]);
                if(d > dmax)
                    dmax = d;
            }
        }

        while(rsize < 2*dmax)
            rsize = 2*rsize;
    }

    void LoadBody(Body* p)
    {
        Cell* q; Cell* c;
        int qind,k;
        double qsize, dist2;

        q = root;
        qind = SubIndex(p,q);
        qsize = rsize;

        while(q->subP[qind] != nullptr)
        {
            if(q->subP[qind]->type == Type::Body)
            {
                c = MakeCell();
                for(int k = 0;k<DIM;++k)
                    c->position[k] = q->position[k] + (p->position[k] < q->position[k] ? -qsize:qsize)/4;
                c->subP[SubIndex(static_cast<Body*>(q->subP[qind]),c)] = q->subP[qind];
                q->subP[qind] = c;
            }
            q = static_cast<Cell*>(q->subP[qind]);
            qind = SubIndex(p,q);
            qsize /= 2;
        }

        q->subP[qind] = p;
    }

    int SubIndex(Body* p, Cell* q)
    {
        int ind, k;

        ind = 0;
        for(k = 0;k<DIM;++k)
        {
            if(q->position[k] <= p->position[k])
                ind += NSUB >> (k+1);
        }

        return ind;
    }

    void PropagateInfo(Cell* p, double psize, int lev, std::function<void(Cell*)> initializer, std::function<void(Cell*, double, Cell*)> propagater_cc, std::function<void(Cell*, double, int)> propagater_cb, std::function<void(Cell*, double)> dataProcessor)
    {
        Node* q;

        tdepth = std::max(tdepth,lev);
        initializer(p);
        for(int i = 0;i<NSUB;++i)
        {
            if((q=p->subP[i]) != nullptr)
            {
                if(q->type == Type::Cell)
                    PropagateInfo(static_cast<Cell*>(q), psize/2, lev+1, initializer, propagater_cc, propagater_cb, dataProcessor);
                
                p->update |= q->update;
                if(q->type == Type::Cell)
                    propagater_cc(p,psize,static_cast<Cell*>(q));
                else
                    propagater_cb(p,psize,static_cast<Body*>(q)->id);
            }
        }

        dataProcessor(p,psize);
    }

    void ThreadTree(Node* p, Node* n)
    {
        int ndesc, i;
        Node* desc[NSUB+1],*_p;
        Cell* _cellP;

        p->next = n;
        if(p->type == Type::Cell)
        {
            _cellP = static_cast<Cell*>(p);
            ndesc = 0;
            for(i = 0;i<NSUB;++i)
            {
                _p = _cellP->subP[i];
                if(_p != nullptr)
                    desc[ndesc++] = _p;
            }            
            _cellP->more = desc[0];
            desc[ndesc] = n;
            for(i = 0;i<ndesc;++i)
                ThreadTree(desc[i], desc[i+1]);
        }
    }

    void WalkTree(Node** aptr, Node** nptr, Cell* cptr, int* bptr, int i, Node* p, double psize, double* pmid)
    {
        Node **np, **ap, *q;
        int actSafe;

        if(p->update)
        {
            np = nptr;
            actSafe = actLen - NSUB;
            for(ap = aptr; ap < nptr; ++ap)
            {
                if((*ap)->type == Type::Cell)
                {
                    if(accept(static_cast<Cell*>(*ap),psize,pmid))
                    {
                        cptr->position = (*ap)->position;
                        cptr->info = static_cast<Cell*>(*ap)->info;
                        ++cptr;
                    }else
                    {
                        if(np-active >= actSafe)
                        {
                            std::cout << "dont't have enough memory";
                            std::exit(-1);
                        }
                        for(q = static_cast<Cell*>(*ap)->more;q != (*ap)->next; q = q->next)
                            *np++ = q;
                    }
                }else
                {
                    if((*ap) != p)
                    {
                        bptr[i] = static_cast<Body*>(*ap) -> id;
                        ++i;
                    }
                }
            }

            if(np != nptr)
                WalkSub(nptr,np,cptr,bptr,i,p,psize,pmid);
            else
            {
                if (p->type != Type::Body)                /* must have found a body   */
                    exit(2);
                DoAction(static_cast<Body*>(p),cptr,bptr,i);
            }

        }
    }

    bool isFarFromTarget(Node* c, double psize, double* pmid)
    {
        double dx, farLen;

        farLen = psize + length;

        for(int dim = 0;dim < DIM;++dim)
        {
            dx = c->position[dim] - pmid[dim];
            if(std::abs(dx) > farLen)
                return true;
        }

        return false;
    }

    void WalkSub(Node** nptr, Node** np, Cell* cptr, int* bptr, int i, Node* p, double psize, double* pmid)
    {
        double poff;
        Node* q;
        int k;
        double nmid[DIM];

        poff = psize/4;
        if(p->type == Type::Cell)
        {
            for(q = static_cast<Cell*>(p)->more; q != p->next;q=q->next)
            {
                for(k = 0;k<DIM;++k)
                    nmid[k] = pmid[k] + (q->position[k] < pmid[k] ? -poff:poff);
                WalkTree(nptr,np,cptr,bptr,i,q,psize/2,nmid);
            }
        }else
        {
            for(k = 0;k<DIM;++k)
                nmid[k] = pmid[k] + (p->position[k] < pmid[k] ? -poff:poff);
            WalkTree(nptr,np,cptr,bptr,i,p,psize/2,nmid);
        }
    }

    void DoAction(Body* p0, Cell* cptr, int* bptr, int i)
    {
       int bodyID = p0->id;

        Cell* p;
        for(p = interact;p < cptr;++p)
            bc_action(bodyID,p);
        for(int k = 0;k<i;++k)
            action(bodyID,bptr[k]);
    }
    
public:
    double rsize = 1;//root size
    int tdepth;//木の高さ

    int num;//現在のAタイプ粒子数
    int num2;//現在のBタイプ粒子数

    int num_max;//Aタイプのマックス粒子数
    int num_max2;//Bタイプのマックス粒子数

    BarnesHutTree(const posType& position, unsigned int number, unsigned char _DIM):num(number),num2(0),num_max(number),num_max2(0),DIM(_DIM),NSUB(1<<_DIM)//about particle labeld "i", you must be able to access d-axis position by "position[i][dim]"
    {
        bodies = Allocate_Array<Body>(number);
        for(int i = 0;i<number;++i)
        {
            bodies[i].type = Type::Body;
            bodies[i].update = true;
            bodies[i].position = position[i];
            bodies[i].id = i;
        }
    }

    BarnesHutTree(const posType& position, unsigned int number, const posType& position2, unsigned int number2, unsigned int _DIM):num(number),num2(number2),num_max(number),num_max2(number2),DIM(_DIM),NSUB(1<<_DIM)//about particle labeld "i", you must be able to access d-axis position by "position[i][dim]"
    {
        bodies = Allocate_Array<Body>(number);
        bodies2 = Allocate_Array<Body>(number2);
        for(int i = 0;i<number;++i)
        {
            bodies[i].type = Type::Body;
            bodies[i].update = true;
            bodies[i].position = position[i];
            bodies[i].id = i;
        }

        for(int i = 0;i<number2;++i)
        {
            bodies2[i].type = Type::Body;
            bodies2[i].update = false;
            bodies2[i].position = position2[i];
            bodies2[i].id = i+number;
        }
    }

    ~BarnesHutTree()
    {
        NewTree();
        Cell* c;
        while(freeCell != nullptr)
        {
            c = static_cast<Cell*>(freeCell);
            freeCell = c->next;
            delete c;
        }
        if(freeCell != nullptr)
            delete freeCell;
    
        delete[] bodies;
        
    }
    
    /**
     * @brief make tree
     * 
     * remember the fact that all of edge of tree is either nullptr or body
     * @param initializer   initialize additional cell info
     * @param propagater_cc the function for when cell B is a child of cell A, cell A gets info from cell B and updates additional cell info
     * @param propagater_cb the function for when body is a child of cell, cell gets info from body and update additional cell info
     * @param dataProcessor the function for when info propagation is over, adjust data of cell
     */
    void MakeTree(std::function<void(Cell*)> initializer, std::function<void(Cell*, double, Cell*)> propagater_cc, std::function<void(Cell*, double, int)> propagater_cb, std::function<void(Cell*, double)> dataProcessor)
    {
        Body* p;

        NewTree();
        root = MakeCell();
        root->ClearPosition(DIM);
        ExpandBox();
        for(p = bodies;p<bodies+num;++p)
            LoadBody(p);
        
        tdepth = 0;
        PropagateInfo(root, rsize, 0, initializer, propagater_cc, propagater_cb, dataProcessor);
        ThreadTree(root,nullptr);
        
    }

    void MakeTree(int number, int number2,std::function<void(Cell*)> initializer, std::function<void(Cell*, double, Cell*)> propagater_cc, std::function<void(Cell*, double, int)> propagater_cb, std::function<void(Cell*, double)> dataProcessor)
    {
        Body* p;

        num = number;
        num2 = number2;

        NewTree();
        root = MakeCell();
        root->ClearPosition(DIM);
        ExpandBox();
        for(p = bodies;p<bodies+num;++p)
            LoadBody(p);

        for(p = bodies2;p<bodies2+num2;++p)
            LoadBody(p);
        
        tdepth = 0;
        PropagateInfo(root, rsize, 0, initializer, propagater_cc, propagater_cb, dataProcessor);
        ThreadTree(root,nullptr);
        
    }
    
    /**
     * @brief 
     * 
     * @param _accept 
     * @param _action 
     * @param BC_action 
     */
    void FindNeighborParticleAndAction(std::function<bool(Cell* c, double psize, double* pmid)> _accept, std::function<void(int,int)> _action, std::function<void(int,Cell*)> BC_action)
    {
        double rmid[DIM];
        actLen = 20*216*tdepth;
        accept = _accept;
        action = _action; //BB_action
        bc_action = BC_action;

        active = Allocate_Array<Node*>(actLen);
        interact = Allocate_Array<Cell>(actLen);

        int interact_body[actLen];

        active[0] = root;
        for(int dim = 0;dim<DIM;++dim)
            rmid[dim] = 0;
        WalkTree(active, active+1, interact, interact_body, 0, root, rsize, rmid);

        delete[] active;
        delete[] interact;
    }
};