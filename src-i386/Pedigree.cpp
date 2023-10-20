#include "Pedigree.h"

iNode::iNode(const int i,const int d,const int s)
{
  ind = i;
  dam = d;
  sire = s;
  f = 0.0;
  iNum = -1;
  getA = false;
  gen = max_gen;
  nOff = 0;
  act = -1; // indicates that the individual is not being visited by a recursive function at this moment
  loop = -1; // indicates that the individual is involved in a loop, should check this individual
}

pedigree::~pedigree()
{
  iterator it;
  iNode *pInd;
  for(it = begin();it!=end();it++)
    {
      pInd = (*it);
      delete pInd;
    }
}

double pedigree::getRij(const int dam,const int sire)
{
  if(dam == 0||sire==0){return 0.0;}
  map<const Mij,double>::iterator cit;
  Mij cij;
  int i = dam,j = sire;
  if(i>j){i = sire;j = dam;}

  cij.i = i;
  cij.j = j;
  int dj,sj;
  const iNode *ind;
  double x;
  cit = C.find(cij);
  if(C.find(cij)!=C.end()){
    x = C[cij];
    return x;
  }
  else if(i==j)
    {
      ind = (*this)[i-1];
      x = 0.5*(1+ind->f);
      C[cij] = x;
      return x;
    }
  else
    {
      ind = (*this)[j-1];
      dj = ind->dam;
      sj = ind->sire;
      x = 0.5*(getRij(i,dj) + getRij(i,sj));
      C[cij] = x;
      return(x);
    }
}

void pedigree::makeA()
{
  iterator it,jt;
  iNode *pInd,*poInd;
  int i,j;
  Mij Aij;
  double x;
  for(it = begin();it!=end();it++)
    {
      pInd = (*it);
      if(pInd->getA)
	{
	  i = pInd->ind;
	  Aij.i = Aij.j = i;
	  A[Aij] = 1 + pInd->f;
	  for(jt = begin();jt!=it;jt++)
	    {
	      poInd = (*jt);
	      if(poInd->getA)
		{
		  j = poInd->ind;
		  Aij.j = j;
		  x = 2*getRij(i,j);
		  A[Aij] = x;
		}
	    }
	}
    }
}

void pedigree::writeA(string aFile)
{
  ofstream ainvFile;
  ainvFile.open(aFile.c_str());
  ainvFile.setf(ios_base::left,ios_base::adjustfield);
  map<const Mij,double>::iterator it;
  Mij Aij;
  for(it = A.begin();it!=A.end();it++)
    {
      Aij = (*it).first;
      ainvFile<<Aij.i<<" "<<Aij.j<<" "<<(*it).second<<endl;
    }
  ainvFile.close();
}

void pedigree::trimPedigree()
{
  iterator it;
  iNode *pInd;
  for(it = begin();it!=end();it++)
    {
      pInd = (*it);
      if(pInd->iNum == 0){codeTrim(pInd);} // betekent dat het dier data heeft we er zijn geweest
    }
  C.clear();
}

void pedigree::codeTrim(iNode *pInd)
{
  iNode *pDam,*pSire;
  int gen;
  gen = pInd->gen + 1; 		// ook hier lopen generaties vooruit
  if(pInd->dam != 0)
    {
      pDam = (*this)[pInd->dam - 1];
      if(pDam->gen > gen){pDam->gen = gen;}
      if(pDam->iNum == -1){codeTrim(pDam);} // betekent nog niet gecodeerd
    }
  if(pInd->sire != 0)
    {
      pSire = (*this)[pInd->sire - 1];
      if(pSire->gen > gen){pSire->gen = gen;}
      if(pSire->iNum == -1){codeTrim(pSire);} // betekent nog niet gecodeerd
    }
  pInd->iNum = 1;
}

void pedigree::countOff()
{
  reverse_iterator rit;
  iNode *pInd,*pDam,*pSire;
  for(rit = rbegin();rit<rend();rit++)
    {
      pInd = (*rit);
      if(pInd->dam != 0){
	pDam = (*this)[pInd->dam - 1];
	pDam->nOff = pDam->nOff + pInd->nOff + 1; //noff of this individual plus itself
      }
      if(pInd->sire != 0){
	pSire = (*this)[pInd->sire - 1];
	pSire->nOff = pSire->nOff + pInd->nOff + 1; //noff of this individual plus itself
      }
    }
}


void pedigree::codePedigree()
{
  COUNT = 0;
  iterator it;
  iNode *pInd;
  for(it = begin();it!=end();it++)
    {
      pInd = (*it);
      if(pInd->iNum == -1){code(pInd);}
    }
  C.clear();
}

void pedigree::code(iNode *pInd)
{
  iNode *pDam,*pSire;
  int gDam = 0,gSire = 0;
  if(pInd->act==-1)
    {
      pInd->act = 1;
      if(pInd->dam != 0)
	{
	  pDam = (*this)[pInd->dam - 1];
	  if(pDam->iNum == -1){code(pDam);}
	  gDam = pDam->gen + 1;
	}

      if(pInd->sire != 0)
	{
	  pSire = (*this)[pInd->sire - 1];
	  if(pSire->iNum == -1){code(pSire);}
	  gSire = pSire->gen + 1;
	}
      COUNT++;
      pInd->iNum = COUNT;
      pInd->gen = gDam;
      if(gSire > gDam){pInd->gen = gSire;}
      //cout<<pInd->gen<<endl;
      pInd->act = -1;
    }
  else
    {
      pInd->loop = 1;
    }
}

// extern function to order a pedigree
extern "C"{
  void orderPed(int *ind,int *dam,int *sire,int *n,int *order)
  {
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
      }
    ped->codePedigree();
    for(i = 0;i<*n;i++)
      {
	if((*ped)[i]->loop==1)
	  {
	    order[i] = -1;
	  }
	else
	  order[i] = (*ped)[i]->iNum;
      }
    delete ped;
  }
}


//maak Ainv matrix
void pedigree::makeAinv()
{
  iterator it;
  const iNode *pInd,*pDam,*pSire;
  int i,d,s;
  Mij Aij;
  double D;
  for(it = begin();it!=end();it++)
    {
      pInd = (*it);
      pDam = pSire = NULL;
      if(pInd->dam > 0){pDam = (*this)[pInd->dam - 1];}
      if(pInd->sire > 0){pSire = (*this)[pInd->sire - 1];}

      if(pDam!=NULL&&pSire!=NULL)
	{
	  D = 4.0/( 2.0 - pDam->f - pSire->f);

	  i = pInd->iNum;d = pDam->iNum;s = pSire->iNum;
	  if(s>d)
	    {
	      d = pSire->iNum;
	      s = pDam->iNum;
	    }
	  Aij.i = Aij.j = i;
	  Ainv[Aij] += D;
	  Aij.i = i;Aij.j = d;
	  Ainv[Aij] += -0.5*D;
	  Aij.i = i;Aij.j = s;
	  Ainv[Aij] += -0.5*D;
	  Aij.i = Aij.j = d;
	  Ainv[Aij] += 0.25*D;
	  Aij.i = Aij.j = s;
	  Ainv[Aij] += 0.25*D;
	  Aij.i = d;Aij.j = s;
	  Ainv[Aij] += 0.25*D;
	}
      else if(pDam==NULL&&pSire!=NULL)
	{
	  D = 4/(3 - pSire->f);
	  i = pInd->iNum;s = pSire->iNum;
	  Aij.i = Aij.j = i;
	  Ainv[Aij] += D;
	  Aij.i = i;Aij.j = s;
	  Ainv[Aij] += -0.5*D;
	  Aij.i = Aij.j = s;
	  Ainv[Aij] += 0.25*D;
	}
      else if(pDam!=NULL&&pSire==NULL)
	{
	  D = 4/(3 - pDam->f);
	  i = pInd->iNum;d = pDam->iNum;
	  Aij.i = Aij.j = i;
	  Ainv[Aij] += D;
	  Aij.i = i;Aij.j = d;
	  Ainv[Aij] += -0.5*D;
	  Aij.i = Aij.j = d;
	  Ainv[Aij] += 0.25*D;
	}
      else
	{
	  i = pInd->iNum;
	  Aij.i = Aij.j = i;
	  Ainv[Aij] += 1;
	}
    }
}

void pedigree::writeAinv(string aFile)
{
  ofstream ainvFile;
  ainvFile.open(aFile.c_str());
  ainvFile.setf(ios_base::left,ios_base::adjustfield);
  map<const Mij,double>::iterator it;
  Mij Aij;
  for(it = Ainv.begin();it!=Ainv.end();it++)
    {
      Aij = (*it).first;
      ainvFile<<Aij.i<<" "<<Aij.j<<" "<<(*it).second<<endl;
    }
  ainvFile.close();
}


// externe functies
extern "C"{
  void countOff(int *ind,int *dam,int *sire,int *n,int *nOff)
  {
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
      }
    ped->codePedigree();
    ped->countOff();
    for(i = 0;i<*n;i++)
      {
	nOff[i] = (*ped)[i]->nOff;
      }
    delete ped;
  }
}

extern "C"{
  void countGen (int *ind,int *dam,int *sire,int *n,int *gen)
  {
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
      }
    ped->codePedigree();
    for(i = 0;i<*n;i++)
      {
	gen[i] = (*ped)[i]->gen;
      }
    delete ped;
  }
}

extern "C"{
  void trimPed(int *ind,int *dam,int *sire,int *data,int *ngenback,int *n)
  {
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
	if(data[i]==1){
	  IND->iNum = 0;
	  IND->gen = 0; // als dier data heeft dan gen = 0
	}
      }
    ped->trimPedigree();
    for(i = 0;i<*n;i++)
      {
	if((*ped)[i]->iNum ==1){
	  if((*ped)[i]->gen<=(*ngenback)){data[i] = 1;}
	} // hoeft geen rekening met NULL te houden, iNum van deze dieren is -1
	else {data[i] = 0;}
      }
    delete ped;
  }
}

extern "C"{
  void getAinv(int *ind,int *dam,int *sire,int *n)
  {
    string aFile = "Ainv.txt";
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
	IND->f = ped->getRij(IND->dam,IND->sire);
	IND->iNum = i+1;
      }
    ped->makeAinv();
    ped->writeAinv(aFile);
    delete ped;
  }
}

extern "C" {
  void getA (int *ind,int *dam,int *sire,int *n,int *which)
  {
    string aFile = "A.txt";
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
	IND->f = ped->getRij(IND->dam,IND->sire);
	IND->iNum = i+1;
	if(which[i] == 1){IND->getA = true;}
      }
    ped->makeA();
    ped->writeA(aFile);
    delete ped;
  }
}

//geeft inteeltcoefficienten weer
extern "C" {
  void calcInbreeding (int *ind,int *dam,int *sire,int *n,double *F)
  {
    pedigree *ped;
    iNode *IND;
    int i;
    ped = new pedigree();
    for(i = 0;i<*n;i++)
      {
	IND = new iNode(ind[i],dam[i],sire[i]);
	ped->push_back(IND);
	IND->f = ped->getRij(IND->dam,IND->sire);
	F[i] = IND->f;
      }
    delete ped;
  }
}







