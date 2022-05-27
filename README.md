Given the coordinates of all particles, the space is recursively divided so that only one particle can fit into the divided subspace. By doing so, you can search for another particle within a radius of r from one particle.

You can use this program as a part of NBody-code, sph and collision detection!

## Example

If you just want to search for another particle within a radius of r from one particle.

```c++
    struct AdditionalInfo
    {
        
    };

    const unsigned int number = 1000; //number of particles
    const unsigned int DIM    = 3;    //dimention of space
    double position[number][DIM];     //position of all particles(x,y,z)

    using Tree = BarnesHutTree<double(*)[DIM],AdditionalInfo>; 
    //double(*)[DIM] is a type of two-dimensional array
    Tree tree(position,number,DIM);

    int main()
    {
        while(true)
        {
            tree.MakeTree([](Tree::Cell* p){},
                          //Functor to initialize Additional info of Cell p.
                                                
                          [](Tree::Cell* p,real psize,Tree::Cell* q){}, 
                          //Functor to allow Cell p to get information aboue Cell q which is a child of Cell p. psize is length of a side of Cell p.
                                                    
                          [](Tree::Cell* p, real psize, int id){}, 
                          //Functor to allow Cell p to get information aboue Body id which is a child of Cell p. psize is length of a side of Cell p.
                                                
                          [](Tree::Cell* p,real psize){}); 
                          //Functor to allow cell p to complete data after get all of data form all of chiled of Cell p. psize is length of a side of Cell p.

            
            double r_criterion = 1.0;
            tree.FindNeighborParticleAndAction([&r_criterion](Tree::Cell* c, real psize, real* pmid) 
                                                {
                                                    
                                                    real dx, farLen;

                                                    farLen = psize + r_criterion;

                                                    for(int dim = 0;dim < DIM;++dim)
                                                    {
                                                        dx = c->position[dim] - pmid[dim];
                                                        if(std::abs(dx) > farLen)
                                                            return true;
                                                    }

                                                    return false;
                                                },
                                                [this](int i, int j)     { std::cout << "body " << j << "is most likely within a radius r from body " << i << "\n";},
                                                
                                                [](int i, Tree::Cell* p) { std::cout << "All region of Cell p is outside a radius r from body i"     ;}
                                                )

            UpdatePosition();
        }
    }
```

If you want to implement NBody code,

```c++
    struct AdditionalInfo
    {
        double mass = 0;   // sum of all particles' mass inside the cell
        double rcrit2 = 0; // cellSize/theta
        double cmpos[3];   // center of mass  
    };

    const unsigned int number = 1000; //number of particles
    const unsigned int DIM    = 3;    //dimention of space
    double position[number][DIM];     //position of all particles(x,y,z)
    double mass[number];

    using Tree = BarnesHutTree<double(*)[DIM],AdditionalInfo>; 
    //double(*)[DIM] is a type of two-dimensional array
    Tree tree(position,number,DIM);

    int main()
    {
        while(true)
        {
            tree.MakeTree([&number,&DIM,&position,&tree](Tree::Cell* p) 
                                            {
                                                p->info.mass = 0;
                                                p->info.rcrit2 = 0;
                                                for(int d = 0;d<DIM;++d)
                                                    p->info.cmpos[d] = 0;
                                            }, 
                        //Functor to initialize Additional info of Cell p.
                        [&number,&DIM,&position,&tree](Tree::Cell* p,real psize,Tree::Cell* q)
                                        {    
                                            p->info.mass += q->info.mass;
                                            for(int d = 0;d<DIM;++d)
                                                p->info.cmpos[d] += q->position[d]*q->info.mass;
                                        },
                        //Functor to allow Cell p to get information aboue Cell q which is a child of Cell p. psize is length of a side of Cell p.
                        [&number,&DIM,&position,&tree](Tree::Cell* p, real psize, int id)
                                        {
                                            p->info.mass += mass[id];
                                            for(int d = 0;d<DIM;++d)
                                                p->info.cmpos[d] += position[id][d]*mass[id];
                                        },
                        //Functor to allow Cell p to get information aboue Body id which is a child of Cell p. psize is length of a side of Cell p.
                        [&number,&DIM,&position,&tree](Tree::Cell* p,real psize) 
                                        {
                                            if(p->info.mass > 0)
                                            {
                                                for(int d = 0;d<DIM;++d)
                                                    p->info.cmpos[d] /= p->info.mass;
                                            }else
                                            {
                                                for(int d = 0;d<DIM;++d)
                                                    p->info.cmpos[d] = p->position[d];
                                            }

                                            if(theta == 0)
                                            {
                                                p->info.rcrit2 = 2*tree.rsize;
                                                p->info.rcrit2 *= p->info.rcrit2;
                                            }else
                                            {
                                                p->info.rcrit2 = psize/theta;
                                                p->info.rcrit2 *= p->info.rcrit2;
                                            }

                                            for(int d = 0;d<DIM;++d)
                                            {
                                                p->position[d] = p->info.cmpos[d];
                                            }

                                        }
                            //Functor to allow cell p to complete data after get all of data form all of chiled of Cell p. psize is length of a side of Cell p.
                            );


                    tree.FindNeighborParticleAndAction(
                        [this](Tree::Cell* c, real psize, real* pmid)
                        {
                            
                            real dmax, dsq, dk;
                            int k;
                        
                            dmax = psize;                               /* init maximum distance    */
                            dsq = 0.0;                                  /* and squared min distance */
                            for (k = 0; k < this->DIM; k++) {                /* loop over space dims     */
                                dk = c->position[k] - pmid[k];               /* form distance to midpnt  */
                                if (dk < 0)                             /* and get absolute value   */
                                    dk = - dk;
                                if (dk > dmax)                          /* keep track of max value  */
                                    dmax = dk;
                                dk -= ((real) 0.5) * psize;             /* allow for size of cell   */
                                if (dk > 0)
                                    dsq += dk * dk;                     /* sum min dist to cell ^2  */
                            }
                            return (dsq > c->info.rcrit2 &&                  /* test angular criterion   */
                                    dmax > ((real) 1.5) * psize);     /* and adjacency criterion  */
                        },

                        [this](int i, int j)
                        {
                            std::cout << "interaction with particle " << i << "and partilce " << j << "\n";
                        },

                        [this](int i, Tree::Cell* p)
                        {
                            std::cout << "interaction with particle " << i << "and cell " << p << "\n";
                            //you can accese p->info.mass,p->info.rcrit2,p->info.cmpos,p->position
                        })

            UpdatePosition();
        }
    }
```

## How to use

This program is header only. You can just use BarnesHutTree by just including.

```c++
    #include "BarnesHutTree.hpp"
```

In BarnesHutTree, we call subspace cell, and particle body.
These are what you need at least.

```c++
    struct AdditionalInfo
    {
        
    };

    const unsigned int number = 1000; //number of particles
    const unsigned int DIM    = 3;    //dimention of space
    double position[number][DIM];     //position of all particles(x,y,z)

    using Tree = BarnesHutTree<double(*)[DIM],AdditionalInfo>; 
    //double(*)[DIM] is a type of two-dimensional array
    Tree tree(position,number,DIM);
```

AdditionalInfo is a struct for Additional information you want to store into a cell.

If you just want to search for another particle within a radius of r from one particle, your AdditionalInfo is going to be like this:

```c++
    struct AdditionalInfo
    {
        
    }; 
```

If you want to implement NBody code, your AdditionalInfo is going to be like this:

```c++
    struct AdditionalInfo
    {
        double mass = 0;   // sum of all particles' mass inside the cell
        double rcrit2 = 0; // cellSize/theta
        double cmpos[3];   // center of mass  
    };
```

I think you already know that you can acces pointer to the position of particle-i by this

```c++
    double* position_i = position[i];
    double x           = position_i[0];
    double y           = position_i[1];
    double z           = position_i[2];
```

If your custom array class also can do same thing, you can use your custome array class as your strage class of all particles.

```c++
    YourCustomArrayClass<double> position(number,DIM);     //position of all particles(x,y,z)

    double* position_i = position[i];
    double x           = position_i[0];
    double y           = position_i[1];
    double z           = position_i[2];

    using Tree = BarnesHutTree<YourCustomArrayClass<double>,AdditionalInfo>; 
    //double(*)[DIM] is a type of two-dimensional array
    Tree tree(position,number,DIM);
```

next thing you have to do is this below

```c++
    tree.MakeTree([](Tree::Cell* p){},//Functor to initialize Additional info of Cell p.
                                             
                  [](Tree::Cell* p,real psize,Tree::Cell* q){}, //Functor to allow Cell p to get information aboue Cell q which is a child of Cell p. psize is length of a side of Cell q.
                                            
                  [](Tree::Cell* p, real psize, int id){}, //Functor to allow Cell p to get information aboue Body id which is a child of Cell p. psize is length of a side of Cell q.
                                        
                  [](Tree::Cell* p,real psize){}); //Functor to allow cell p to complete data after get all of data form all of chiled of Cell p.
```

If you just want to search for another particle within a radius of r from one particle, you dont have to write anything in functor above.

If you want to implement NBody code, you have to write functor like this below.
Note that you can acces position of center of Cell p by p->position.

```c++
    tree.MakeTree([](Tree::Cell* p) 
                                    {
                                        p->info.mass = 0;
                                        p->info.rcrit2 = 0;
                                        for(int d = 0;d<DIM;++d)
                                            p->info.cmpos[d] = 0;
                                    }, 
                  [](Tree::Cell* p,real psize,Tree::Cell* q)
                                {    
                                    p->info.mass += q->info.mass;
                                    for(int d = 0;d<DIM;++d)
                                        p->info.cmpos[d] += q->position[d]*q->info.mass;
                                },
                  [](Tree::Cell* p, real psize, int id)
                                {
                                    p->info.mass += mass[id];
                                    for(int d = 0;d<DIM;++d)
                                        p->info.cmpos[d] += position[id][d]*mass[id];
                                },
                  [](Tree::Cell* p,real psize) 
                                {
                                    if(p->info.mass > 0)
                                    {
                                        for(int d = 0;d<DIM;++d)
                                            p->info.cmpos[d] /= p->info.mass;
                                    }else
                                    {
                                        for(int d = 0;d<DIM;++d)
                                            p->info.cmpos[d] = p->position[d];
                                    }

                                    if(theta == 0)
                                    {
                                        p->info.rcrit2 = 2*tree.rsize;
                                        p->info.rcrit2 *= p->info.rcrit2;
                                    }else
                                    {
                                        p->info.rcrit2 = psize/theta;
                                        p->info.rcrit2 *= p->info.rcrit2;
                                    }

                                    for(int d = 0;d<DIM;++d)
                                    {
                                        p->position[d] = p->info.cmpos[d];
                                    }

                                });
```

next thing you have to do is this below

```c++
    tree.FindNeighborParticleAndAction([](Tree::Cell* c, real psize, real* pmid) 
                                                {
                                                    if(a)
                                                        return false;
                                                    if(b)
                                                        return true;
                                                },
                                                //Functor to check if Cell c pass the criterion. pmid, psize is another cell of position.
                                                
                                                [this](int i, int j)     { },
                                                
                                                [](int i, Tree::Cell* p) { }
                                                )
```

If you just want to search for another particle within a radius of r from one particle

```c++
    double r_criterion = 1.0;
    tree.FindNeighborParticleAndAction([&r_criterion](Tree::Cell* c, real psize, real* pmid) 
                                        {
                                            
                                            real dx, farLen;

                                            farLen = psize + r_criterion;

                                            for(int dim = 0;dim < DIM;++dim)
                                            {
                                                dx = c->position[dim] - pmid[dim];
                                                if(std::abs(dx) > farLen)
                                                    return true;
                                            }

                                            return false;
                                        },
                                        [this](int i, int j)     { std::cout << "body " << j << "is within a radius r from body " << i << "\n";},
                                        
                                        [](int i, Tree::Cell* p) { std::cout << "All region of Cell p is outside a radius r from body i"     ;}
                                        )
```

If you want to implement NBody code, you have to write functor like this below.

```c++
    tree.FindNeighborParticleAndAction(
                        [this](Tree::Cell* c, real psize, real* pmid)
                        {
                            
                            real dmax, dsq, dk;
                            int k;
                        
                            dmax = psize;                               /* init maximum distance    */
                            dsq = 0.0;                                  /* and squared min distance */
                            for (k = 0; k < this->DIM; k++) {                /* loop over space dims     */
                                dk = c->position[k] - pmid[k];               /* form distance to midpnt  */
                                if (dk < 0)                             /* and get absolute value   */
                                    dk = - dk;
                                if (dk > dmax)                          /* keep track of max value  */
                                    dmax = dk;
                                dk -= ((real) 0.5) * psize;             /* allow for size of cell   */
                                if (dk > 0)
                                    dsq += dk * dk;                     /* sum min dist to cell ^2  */
                            }
                            return (dsq > c->info.rcrit2 &&                  /* test angular criterion   */
                                    dmax > ((real) 1.5) * psize);     /* and adjacency criterion  */
                        },

                        [this](int i, int j)
                        {
                            std::cout << "interaction with particle " << i << "and partilce " << j << "\n";
                        },

                        [this](int i, Tree::Cell* p)
                        {
                            std::cout << "interaction with particle " << i << "and cell " << p << "\n";
                            //you can accese p->info.mass,p->info.rcrit2,p->info.cmpos,p->position
                        })
```
