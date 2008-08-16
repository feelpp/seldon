#ifndef FILE_SELDON_FUNCTIONS_ARRAYS_CXX

namespace Seldon
{
  
  
  ////////////
  //  SORT  //
  
  
  //! Intermediary function used for quick sort algorithm.
  template<class T, class Storage, class Allocator>
  int PartitionQuickSort(int m, int n,
			 Vector<T, Storage, Allocator>& t)
  {
    T temp, v;
    v = t(m);
    int i = m - 1;
    int j = n + 1;
    
    while (true)
      {
	do
	  {
	    j--;
	  }
	while (t(j) > v);
	
	do
	  {
	    i++;
	  }
	while (t(i) < v);
	
	if (i < j)
	  {
	    temp = t(i);
	    t(i) = t(j);
	    t(j) = temp;
	  }
	else
	  {
	    return j;
	  }
      }
  }
  
  
  //! Vector 't' is sorted by using QuickSort algorithm.
  /*!
    Sorts array 't' between position m and n.
  */
  template<class T, class Storage, class Allocator>
  void QuickSort(int m, int n,
		 Vector<T, Storage, Allocator>& t)
  {
    if (m < n)
      {
	int p = PartitionQuickSort(m, n, t);
	QuickSort(m, p, t);
	QuickSort(p+1, n, t);
      }
  }
  
  
  //! Intermediary function used for quick sort algorithm.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  int PartitionQuickSort(int m, int n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2)
  {
    T1 temp1, v;
    T2 temp2;
    v = t1(m);
    int i = m - 1;
    int j = n + 1;
    
    while (true)
      {
	do
	  {
	    j--;
	  }
	while (t1(j) > v);
	do
	  {
	    i++;
	  }
	while (t1(i) < v);
	
	if (i < j)
	  {
	    temp1 = t1(i);
	    t1(i) = t1(j);
	    t1(j) = temp1;
	    temp2 = t2(i);
	    t2(i) = t2(j);
	    t2(j) = temp2;
	  }
	else
	  {
	    return j;
	  }
      }
  }
  
  
  //! Vector 't1' is sorted by using QuickSort algorithm.
  /*!  Sorts array 't2' between position m and n, the sorting operation
    affects vector 't2'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void QuickSort(int m, int n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2)
  {
    if (m < n)
      {
	int p = PartitionQuickSort(m, n, t1, t2);
	QuickSort(m, p, t1, t2);
	QuickSort(p+1, n, t1, t2);
      }
  }
  
  
  //! Intermediary function used for quick sort algorithm.
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  int PartitionQuickSort(int m, int n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2,
			 Vector<T3, Storage3, Allocator3>& t3)
  {
    T1 temp1, v;
    T2 temp2;
    T3 temp3;
    v = t1(m);
    int i = m - 1;
    int j = n + 1;
    
    while (true)
      {
	do
	  {
	    j--;
	  }
	while (t1(j) > v);
	
	do
	  {
	    i++;
	  }
	while (t1(i) < v);
	
	if (i < j)
	  {
	    temp1 = t1(i);
	    t1(i) = t1(j);
	    t1(j) = temp1;
	    temp2 = t2(i);
	    t2(i) = t2(j);
	    t2(j) = temp2;
	    temp3 = t3(i);
	    t3(i) = t3(j);
	    t3(j) = temp3;
	  }
	else
	  {
	    return j;
	  }
      }
  }
  
  
  //! Vector 't1' is sorted by using QuickSort algorithm.
  /*!  Sorts array 't1' between position m and n, the sorting operation
    affects vectors 't2' and 't3'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void QuickSort(int m, int n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2,
		 Vector<T3, Storage3, Allocator3>& t3)
  {
    if (m < n)
      {
	int p = PartitionQuickSort(m, n, t1, t2, t3);
	QuickSort(m, p, t1, t2, t3);
	QuickSort(p+1, n, t1, t2, t3);
      }
  }
  
  
  //! Vector 'tab1' is sorted by using MergeSort algorithm.
  /*!
    Sorts array 'tab1' between position m and n.
  */
  template<class T, class Storage, class Allocator>
  void MergeSort(int m, int n, Vector<T, Storage, Allocator>& tab1)
  {
    if (m <= n)
      return;
    
    int inc = 1, ind = 0, current, i, j, sup;
    Vector<T, Storage, Allocator> tab1t(m-n+1);
    // We perform a merge sort with a reccurence
    // inc = 1, 2, 4, 8, ...
    while (inc < n)
      {
	// i is the first index of the sub array of size 2*inc.
	// We make a loop on these sub arrays.
	// Each sub array is divided in two sub arrays of size inc.
	// We merge these two sub arrays in one.
	for (i = m; i < n - inc; i += 2 * inc)
	  {
	    ind = i;
	    current = i + inc; // Index of the second sub array.
	    sup = i + 2 * inc; // End of the merged array.
	    if (sup >= n)
	      sup = n;
	    j = i;
	    // We make a loop on values of the first sub array.
	    while (j < i + inc)
	      {
		// If the second sub array has still elements not sorted.
		if (current < sup)
		  {
		    // We insert element of the second sub array in the
		    // merged array until tab1(j) < tab1(current).
		    while (current < sup && tab1(j) > tab1(current))
		      {
			tab1t(ind-m) = tab1(current);
			current++;
			ind++;
		      }
		    // We insert the element of the first sub array now.
		    tab1t(ind-m) = tab1(j);
		    ind++;
		    j++;
		    // If first sub array is sorted, we insert all remaining
		    // elements of the second sub array.
		    if (j == i + inc)
		      {
			for (j = current; j < sup; j++)
			  {
			    tab1t(ind-m) = tab1(j);
			    ind++;
			  }
		      }
		  }
		else
		  {
		    // If the second sub array is sorted, we insert all
		    // remaining elements of the first sub array.
		    for (current = j; current < i + inc; current++)
		      {
			tab1t(ind-m) = tab1(current);
			ind++;
		      }
		    j = current+1;
		  }
	      }
	  }
	
	for (i = m; i < ind; i++)
	  tab1(i) = tab1t(i-m);
	
	inc = 2*inc;
      }
  }
  
  
  //! Vector 'tab1' is sorted by using MergeSort algorithm.
  /*!  Sorts array 'tab1' between position m and n. The sort operation affects
    'tab2'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void MergeSort(int m, int n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2)
  {
    if (m <= n)
      return;
    
    int inc = 1, ind = 0, current, i, j, sup;
    Vector<T1, Storage1, Allocator1> tab1t(m-n+1);
    Vector<T2, Storage2, Allocator2> tab2t(m-n+1);
    
    while (inc < n)
      {
	for (i = 0; i < n - inc; i += 2 * inc)
	  {
	    ind = i;
	    current = i + inc;
	    sup = i + 2 * inc;
	    if (sup >= n)
	      sup = n;
	    j = i;
	    while (j < i + inc)
	      {
		if (current < sup)
		  {
		    while (current < sup && tab1(j) > tab1(current))
		      {
			tab1t(ind-m) = tab1(current);
			tab2t(ind-m) = tab2(current);
			current++;
			ind++;
		      }
		    tab1t(ind-m) = tab1(j);
		    tab2t(ind-m) = tab2(j);
		    ind++;
		    j++;
		    if (j == i + inc)
		      {
			for (j = current; j < sup; j++)
			  {
			    tab1t(ind-m) = tab1(j);
			    tab2t(ind-m) = tab2(j);
			    ind++;
			  }
		      }
		  }
		else
		  {
		    for (current = j; current < i + inc; current++)
		      {
			tab1t(ind-m) = tab1(current);
			tab2t(ind-m) = tab2(current);
			ind++;
		      }
		    j = current + 1;
		  }
	      }
	  }
	for (i = 0; i < ind; i++)
	  {
	    tab1(i) = tab1t(i-m);
	    tab2(i) = tab2t(i-m);
	  }
	inc = 2 * inc;
      }
  }
  
  
  //! Vector 'tab1' is sorted by using MergeSort algorithm.
  /*!  Sorts array 'tab1' between position m and n. The sort operation affects
    'tab2' and 'tab3'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void MergeSort(int m, int n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2,
		 Vector<T3, Storage3, Allocator3>& tab3)
  {
    if (m <= n)
      return;
    
    int inc = 1, ind = 0, current, i, j, sup;
    Vector<T1, Storage1, Allocator1> tab1t(n-m+1);
    Vector<T2, Storage2, Allocator2> tab2t(n-m+1);
    Vector<T3, Storage3, Allocator3> tab3t(n-m+1);
    
    while (inc < n)
      {
	for (i = 0; i < n - inc; i += 2 * inc)
	  {
	    ind = i;
	    current = i + inc;
	    sup = i + 2 * inc;
	    if (sup >= n)
	      sup = n;
	    j = i;
	    while (j < i + inc)
	      {
		if (current < sup)
		  {
		    while (current < sup && tab1(j) > tab1(current))
		      {
			tab1t(ind-m) = tab1(current);
			tab2t(ind-m) = tab2(current);
			tab3t(ind-m) = tab3(current);
			current++;
			ind++;
		      }
		    tab1t(ind-m) = tab1(j);
		    tab2t(ind-m) = tab2(j);
		    tab3t(ind-m) = tab3(j);
		    ind++;
		    j++;
		    if (j == i + inc)
		      {
			for (j = current; j < sup; j++)
			  {
			    tab1t(ind-m) = tab1(j);
			    tab2t(ind-m) = tab2(j);
			    tab3t(ind-m) = tab3(j);
			    ind++;
			  }
		      }
		  }
		else
		  {
		    for (current = j; current < i + inc; current++)
		      {
			tab1t(ind-m) = tab1(current);
			tab2t(ind-m) = tab2(current);
			tab3t(ind-m) = tab3(current);
			ind++;
		      }
		    j = current+1;
		  }
	      }
	  }
	for (i = 0; i < ind; i++)
	  {
	    tab1(i) = tab1t(i-m);
	    tab2(i) = tab2t(i-m);
	    tab3(i) = tab3t(i-m);
	  }
	inc = 2 * inc;
      }
  }
  
  
  //! Assembles a sparse vector.
  /*!
    Node are the row numbers of the vector.
    Vect are the values of the vector.
    The function sorts the row numbers and
    adds values corresponding to the same row number.
    For example, if Node = [3 2 2 0], Vect = [1.0 0.4 0.4 -0.3]
    the function will return
    n = 3, Node = [0 2 3], Vect = [-0.3 0.8 1.0].
    \param[inout] n on input, number of row numbers to assemble,
    on output, number of row numbers after assembling.
    \param[inout] Node row numbers.
    \param[inout] Vect value numbers.
    \warning Vectors Node and Vect are not resized.
  */
  template<class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2 >
  void Assemble(int& n, Vector<int, Storage1, Allocator1>& Node,
		Vector<T2, Storage2, Allocator2>& Vect)
  {
    if (n <= 1)
      return;
    
    Sort(n, Node, Vect);
    int prec = Node(0);
    int nb = 0;
    for (int i = 1; i < n; i++)
      if (Node(i) == prec)
	{
	  Vect(nb) += Vect(i);
	}
      else
	{
	  nb++;
	  Node(nb) = Node(i);
	  Vect(nb) = Vect(i);
	  prec = Node(nb);
	}
    n = nb + 1;
  }
  
  
  //! Sorts and removes duplicate entries of a vector.
  /*!
    \param[inout] n on input, number of elements to assemble
    on output, number of elements after assembling.
    \param[inout] Node vector to assemble.
    \warning The vector is not resized.
  */
  template<class T, class Storage1, class Allocator1>
  void Assemble(int& n, Vector<T, Storage1, Allocator1>& Node)
  {
    if (n <= 1)
      return;
    
    Sort(n, Node);
    T prec = Node(0);
    int nb = 1;
    for (int i = 1; i < n; i++)
      if (Node(i)!=prec)
	{
	  Node(nb) = Node(i);
	  prec = Node(nb);
	  nb++;
	}
    n = nb;
  }
  

  //! Sorts and removes duplicate entries of a vector.
  template<class T, class Storage1, class Allocator1>
  void Assemble(Vector<T, Storage1, Allocator1>& Node)
  {
    int nb = Node.GetM();
    Assemble(nb, Node);
    Node.Resize(nb);
  }
  
  
  //! Sorts and removes duplicate entries of a vector.
  /*!
    Sorting operations of 'Node' affects 'Node2'.
  */
  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2)
  {
    if (n <= 1)
      return;
    
    Sort(n, Node, Node2);
    T prec = Node(0);
    int nb = 1;
    for (int i = 1; i < n; i++)
      if (Node(i) != prec)
	{
	  Node(nb) = Node(i);
	  Node2(nb) = Node2(i);
	  prec = Node(nb);
	  nb++;
	}
    n = nb;
  }
  

  //! Sorts 'a' and 'b'.
  template<class T>
  inline void Sort(T& a, T& b)
  {
    if (b < a)
      {
	T temp = a;
	a = b;
	b = temp;
      }
  }
  
  
  //! Sorts 'a', 'b' and 'c'.
  template<class T>
  inline void Sort(T& a, T& b, T& c)
  {
    if (b < a)
      {
	T temp = a;
	a = b;
	b = temp;
      }
    if (c < a)
      {
	T temp = a;
	a = c;
	T temp2 = b;
	b = temp;
	c = temp2;
      }
    else if (c < b)
      {
	T temp = c;
	c = b;
	b = temp;
      }
  }
  
  
  //! Sorts 'i', 'j', 'k' and 'l'.
  template<class T>
  inline void Sort(T& i, T& j, T& k, T& l)
  {
    T i0 = i, j0 = j, k0 = k, l0 = l;
    if (i > j)
      {
	i0 = j;
	j0 = i;
      }
    
    if (k > l)
      {
	k0 = l;
	l0 = k;
      }
    
    if (i0 < k0)
      {
	i = i0;
	if (j0 < l0)
	  {
	    l = l0;
	    if (j0 < k0)
	      {
		j = j0;
		k = k0;
	      }
	    else
	      {
		j = k0;
		k = j0;
	      }
	  }
	else
	  {
	    j = k0;
	    k = l0;
	    l = j0;
	  }
      }
    else
      {
	i = k0;
	if (l0 < j0)
	  {
	    l = j0;
	    if (l0 < i0)
	      {
		j = l0;
		k = i0;
	      }
	    else
	      {
		j = i0;
		k = l0;
	      }
	  }
	else
	  {
	    j = i0;
	    k = j0;
	    l = l0;
	  }
      }
  }
  
  
  //! Sorts vector 'V' between a start position and an end position.
  template<class T, class Storage, class Allocator>
  void Sort(int m, int n, Vector<T, Storage, Allocator>& V)
  {
    QuickSort(m, n, V);
  }
  
  
  //! Sorts vector 'V' between a start position and an end position.
  /*!
    The sorting operation of 'V' affects 'V2'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(int m, int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2)
  {
    QuickSort(m, n, V, V2);
  }
  
  
  //! Sorts vector 'V' between a start position and an end position.
  /*!
    The sorting operation of 'V' affects 'V2' and 'V3'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(int m, int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3)
  {
    QuickSort(m, n, V, V2, V3);
  }
  
  
  //! Sorts the 'n' first elements of 'V'.
  template<class T, class Storage, class Allocator>
  void Sort(int n, Vector<T, Storage, Allocator>& V)
  {
    Sort(0, n-1, V);
  }
  
  
  //! Sorts the 'n' first elements of 'V'.
  /*!
    The sorting operation of 'V' affects 'V2'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2)
  {
    Sort(0, n-1, V, V2);
  }
  
  
  //! Sorts the 'n' first elements of 'V'.
  /*!
    The sorting operation of 'V' affects 'V2' and 'V3'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3)
  {
    Sort(0, n-1, V, V2, V3);
  }
  
  
  //! Sorts vector 'V'.
  template<class T, class Storage, class Allocator>
  void Sort(Vector<T, Storage, Allocator>& V)
  {
    Sort(0, V.GetM()-1, V);
  }
  
  
  //! Sorts vector 'V'.
  /*!
    The sorting operation of 'V' affects 'V2'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2)
  {
    Sort(0, V.GetM()-1, V, V2);
  }
  
  
  //! Sorts vector 'V'.
  /*!
    The sorting operation of 'V' affects 'V2' and 'V3'.
  */
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3)
  {
    Sort(0, V.GetM()-1, V, V2, V3);
  }
  
  
  //  SORT  //
  ////////////
  
  
  //! Appends two vectors X & Y -> X.
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2 >
  void Append(Vector<T, Storage1, Allocator1>& X, int n,
	      const Vector<T, Storage2, Allocator2>& Y)
  {
    int m = X.GetM();
    if (n <= 0)
      return;
    
    int length = m + n;
    Vector<T, Storage1, Allocator1> X_new(length);
    for (int i = 0; i < m; i++)
      X_new(i) = X(i);
      
    for (int i = 0; i < n; i++)
      X_new(m+i) = Y(i);
    
    X.SetData(length, X_new.GetData());
    X_new.Nullify();
  }
  
  
  //! Appends two vectors X & Y -> X.
  template<class T, class Storage1, class Allocator1,
	   class Storage2, class Allocator2 >
  void Append(Vector<T, Storage1, Allocator1>& X,
	      const Vector<T, Storage2, Allocator2>& Y)
  {
    Append(X, Y.GetM(), Y);
  }

  
} // end namespace

#define FILE_SELDON_FUNCTIONS_ARRAYS_CXX
#endif
