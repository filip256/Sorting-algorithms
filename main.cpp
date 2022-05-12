#include <iostream>
#include <vector>
#include <SFML/Graphics.hpp>

#define WINDOW_WIDTH 1800
#define WINDOW_HEIGHT 1000

class Plot
{
	std::vector<sf::Vector2f> m_Points;
	sf::Vector2f m_Offset, m_Scale;
	sf::Color m_Color;
	sf::Font m_Font;
	sf::Text m_Name;

public:
	Plot()
	{
		m_Offset = sf::Vector2f{ 0.0, 0.0 };
		m_Color = sf::Color::Black;
		m_Scale = sf::Vector2f(1.0, 1.0);
	}
	Plot(const sf::Vector2f& offset, const sf::Color& color, const sf::Vector2f& scale, const std::string& name, const sf::Vector2f& namePos)
	{
		m_Offset = offset;
		m_Color = color;
		m_Scale = scale;

		m_Font.loadFromFile("dosfont.ttf");
		m_Name.setFont(m_Font);
		m_Name.setFillColor(m_Color);
		m_Name.setPosition(namePos);
		m_Name.setScale(sf::Vector2f(0.75, 0.75));
		m_Name.setOutlineThickness(1.0);
		m_Name.setOutlineColor(sf::Color::Black);
		m_Name.setString(name);
	}

	void setOffset(const sf::Vector2f& offset) { m_Offset = offset; }
	void setScale(const sf::Vector2f& scale) { m_Scale = scale; }
	void setColor(const sf::Color& color) { m_Color = color; }

	void addPoint(const sf::Vector2f& point) { m_Points.push_back(point); }

	void draw(sf::RenderWindow &window)
	{
		const size_t size = m_Points.size();
		if (size > 1)
		{
			for (int i = 0; i < size - 1; ++i)
			{
				sf::VertexArray line(sf::Lines, 2);
				line[0] = sf::Vertex(sf::Vector2f(m_Points[i].x * m_Scale.x + m_Offset.x, WINDOW_HEIGHT - (m_Points[i].y * m_Scale.y + m_Offset.y)), m_Color);
				line[1] = sf::Vertex(sf::Vector2f(m_Points[i + 1].x * m_Scale.x + m_Offset.x, WINDOW_HEIGHT - (m_Points[i + 1].y * m_Scale.y + m_Offset.y)), m_Color);
				window.draw(line);

				sf::RectangleShape p(sf::Vector2f(3.0, 3.0));
				p.setFillColor(m_Color);
				p.setPosition(sf::Vector2f(line[1].position.x - 1, line[1].position.y - 1));
				window.draw(p);
			}
		}
		m_Name.setFont(m_Font);
		window.draw(m_Name);
	}

};

class Grid
{
	int m_MaxValX, m_MaxValY, m_MinValX, m_MinValY;
	sf::Vector2f m_Size, m_Position;
	sf::VertexArray m_XLine{ sf::Lines, 2 }, m_YLine{ sf::Lines, 2 };
	sf::Font m_Font;

public:
	Grid()
	{
		m_MaxValX = 0;
		m_MaxValY = 0;
		m_Position = sf::Vector2f(0, 0);
		m_Size = sf::Vector2f(0, 0);
	}
	Grid(const sf::Vector2f &position, const sf::Vector2f &size, const int maxX, const int maxY, const int minX, const int minY)
	{
		m_MaxValX = maxX;
		m_MaxValY = maxY;
		m_MinValX = minX;
		m_MinValY = minY;
		m_Position = position;
		m_Size = size;
		m_XLine[0] = sf::Vertex(sf::Vector2f(position.x, WINDOW_HEIGHT - position.y), sf::Color::Black);
		m_YLine[0] = sf::Vertex(sf::Vector2f(position.x, WINDOW_HEIGHT - position.y), sf::Color::Black);
		m_XLine[1] = sf::Vertex(sf::Vector2f(position.x + size.x, WINDOW_HEIGHT - position.y), sf::Color::Black);
		m_YLine[1] = sf::Vertex(sf::Vector2f(position.x, WINDOW_HEIGHT - (position.y + size.y)), sf::Color::Black);

		m_Font.loadFromFile("dosfont.ttf");
	}

	sf::Vector2f getPosition() { return m_Position; }
	sf::Vector2f getSize() { return m_Size; }
	int getMaxX() { return m_MaxValX; }
	int getMaxY() { return m_MaxValY; }
	int getMinX() { return m_MinValX; }
	int getMinY() { return m_MinValY; }

	void draw(sf::RenderWindow &window)
	{
		window.draw(m_XLine);
		window.draw(m_YLine);

		sf::Text text;
		text.setFont(m_Font);
		text.setFillColor(sf::Color::Black);
		text.setScale(0.5, 0.5);

		text.setString("0");
		text.setPosition(sf::Vector2f(m_Position.x - 13.0, WINDOW_HEIGHT - m_Position.y));
		window.draw(text);

		sf::RectangleShape axis(sf::Vector2f(1.0, m_Size.y));
		axis.setFillColor(sf::Color(220, 220, 220, 255));
		for (int i = m_MaxValX / 10; i <= m_MaxValX; i += m_MaxValX / 10)
		{
			text.setString(std::to_string(i + m_MinValX));
			text.setPosition(sf::Vector2f(i * m_Size.x / m_MaxValX + m_Position.x - text.getGlobalBounds().width / 2.0, WINDOW_HEIGHT - m_Position.y + 2.0));
			window.draw(text);

			axis.setPosition(sf::Vector2f(i * m_Size.x / m_MaxValX + m_Position.x, WINDOW_HEIGHT - m_Position.y - m_Size.y));
			window.draw(axis);
		}
		axis.setSize(sf::Vector2f(m_Size.x, 1.0));
		for (int i = m_MaxValY / 10; i <= m_MaxValY; i += m_MaxValY / 10)
		{
			text.setString(std::to_string(i + m_MinValX));
			text.setPosition(sf::Vector2f(m_Position.x - text.getGlobalBounds().width - 10.0, WINDOW_HEIGHT - (i * m_Size.y / m_MaxValY + m_Position.y + text.getGlobalBounds().height)));
			window.draw(text);

			axis.setPosition(sf::Vector2f(m_Position.x, WINDOW_HEIGHT - (i * m_Size.y / m_MaxValY + m_Position.y)));
			window.draw(axis);
		}
	}
};

class Graph
{
	std::vector<Plot> m_Plots;
	Grid m_Grid;
	sf::Font m_Font;
	sf::Text m_Text;

public:
	Graph(const sf::Vector2f &position, const sf::Vector2f &size, const int maxX, const int maxY, const int minX, const int minY)
	{
		m_Grid = Grid(position, size, maxX, maxY, minX, minY);
		m_Font.loadFromFile("dosfont.ttf");
		m_Text.setFont(m_Font);
		m_Text.setFillColor(sf::Color::Black);
		m_Text.setScale(sf::Vector2f(0.75, 0.75));
	}

	void addPointToPlot(const sf::Vector2f& point, const int index)
	{
		if (index < m_Plots.size())
			m_Plots[index].addPoint(point);
	}
	void addPlot(const std::string name, const sf::Color& color)
	{
		const sf::Vector2f size = m_Grid.getSize(), position = m_Grid.getPosition();
		m_Plots.push_back(Plot(position, color, sf::Vector2f(size.x / (m_Grid.getMaxX() - m_Grid.getMinX()), size.y / (m_Grid.getMaxY() - m_Grid.getMinY())), name, sf::Vector2f(position.x + size.x + 20.0, m_Plots.size() * 30 + position.y)));
	}

	void draw(sf::RenderWindow &window)
	{
		m_Text.setFont(m_Font);
		const sf::Vector2f gSize = m_Grid.getSize(), gPosition = m_Grid.getPosition();
		m_Text.setString("Array Size");
		m_Text.setPosition(sf::Vector2f(gPosition.x + gSize.x / 2.0 - m_Text.getGlobalBounds().width / 2.0, WINDOW_HEIGHT - gPosition.y + 20.0));
		window.draw(m_Text);

		m_Text.setString("M\ni\nc\nr\no\ns\ne\nc\no\nn\nd\ns");
		m_Text.setLineSpacing(0.75);
		m_Text.setPosition(sf::Vector2f(gPosition.x - 70.0, WINDOW_HEIGHT - (gPosition.y + gSize.y / 2.0 + m_Text.getGlobalBounds().height / 2.0)));
		window.draw(m_Text);
		m_Text.setLineSpacing(1);

		m_Grid.draw(window);
		const size_t size = m_Plots.size();
		for (size_t i = 0; i < size; ++i)
			m_Plots[i].draw(window);
	}
};

void genRandArray(int arr[], const int n)
{
	for (int i = 0; i < n; ++i)
		//arr[i] = rand() & 65535;
		arr[i] = rand() & 262143;
	std::cout << "rand array generated (" << n << ")\n";
}
void genSortedArray(int arr[], const int n)
{
	for (int i = 0; i < n; ++i)
		arr[i] = i;
	//std::swap(arr[n-6], arr[n-1]);
	//std::swap(arr[0], arr[n - 1]);
	//std::swap(arr[n - 4], arr[n - 2]);
	std::cout << "sorted array generated (" << n << ")\n";
}
void genReversedArray(int arr[], const int n)
{
	for (int i = 0; i < n; ++i)
		arr[i] = n - i;
	std::cout << "reversed array generated (" << n << ")\n";
}

void copyArray(int dest[], int src[], int n)
{
	for (int i = 0; i < n; ++i)
		dest[i] = src[i];
}
void printArray(int arr[], int size)
{
	int i;
	for (i = 0; i < size; i++)
		std::cout << arr[i] << ' ';
	std::cout << '\n';
}

void swap(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}
void selectionSort(int arr[], int n)
{
	int i, j, min_idx;

	// One by one move boundary of unsorted subarray 
	for (i = 0; i < n - 1; i++)
	{
		// Find the minimum element in unsorted array 
		min_idx = i;
		for (j = i + 1; j < n; j++)
			if (arr[j] < arr[min_idx])
				min_idx = j;

		// Swap the found minimum element with the first element 
		swap(&arr[min_idx], &arr[i]);
	}
}
void insertionSort(int arr[], int n)
{
	int i, key, j;
	for (i = 1; i < n; i++)
	{
		key = arr[i];
		j = i - 1;

		while (j >= 0 && arr[j] > key)
		{
			arr[j + 1] = arr[j];
			j = j - 1;
		}
		arr[j + 1] = key;
	}
}
void bubbleSort(int arr[], int n)
{
	int i, j;
	for (i = 0; i < n - 1; i++)

		// Last i elements are already in place
		for (j = 0; j < n - i - 1; j++)
			if (arr[j] > arr[j + 1])
				swap(&arr[j], &arr[j + 1]);
}
void shellSort(int arr[], int n)
{
	for (int gap = n / 2; gap > 0; gap /= 2)
	{
		for (int i = gap; i < n; i += 1)
		{
			int temp = arr[i];
			int j;
			for (j = i; j >= gap && arr[j - gap] > temp; j -= gap)
				arr[j] = arr[j - gap];
			arr[j] = temp;
		}
	}
}
void cocktailSort(int arr[], int n)
{
	bool swapped = true;
	int start = 0;
	int end = n - 1;

	while (swapped)
	{
		swapped = false;
		for (int i = start; i < end; ++i)
		{
			if (arr[i] > arr[i + 1]) {
				std::swap(arr[i], arr[i + 1]);
				swapped = true;
			}
		}

		if (!swapped)
			break;

		swapped = false;
		--end;

		for (int i = end - 1; i >= start; --i)
		{
			if (arr[i] > arr[i + 1]) {
				std::swap(arr[i], arr[i + 1]);
				swapped = true;
			}
		}
		++start;
	}
}
void oddEvenSort(int arr[], int n)
{
	bool isSorted = false; // Initially array is unsorted

	while (!isSorted)
	{
		isSorted = true;

		// Perform Bubble sort on odd indexed element
		for (int i = 1; i <= n - 2; i = i + 2)
		{
			if (arr[i] > arr[i + 1])
			{
				std::swap(arr[i], arr[i + 1]);
				isSorted = false;
			}
		}

		// Perform Bubble sort on even indexed element
		for (int i = 0; i <= n - 2; i = i + 2)
		{
			if (arr[i] > arr[i + 1])
			{
				std::swap(arr[i], arr[i + 1]);
				isSorted = false;
			}
		}
	}

	return;
}
void gnomeSort(int arr[], int n)
{
	int index = 1;

	while (index < n) {
		if (arr[index] >= arr[index - 1])
			index++;
		else {
			std::swap(arr[index], arr[index - 1]);
			index--;
		}
	}
	return;
}

int getNextGap(int gap)
{
	// Shrink gap by Shrink factor
	gap = (gap * 10) / 13;

	if (gap < 1)
		return 1;
	return gap;
}
void combSort(int a[], int n)
{
	// Initialize gap
	int gap = n;

	// Initialize swapped as true to make sure that
	// loop runs
	bool swapped = true;

	// Keep running while gap is more than 1 and last
	// iteration caused a swap
	while (gap != 1 || swapped == true)
	{
		// Find next gap
		gap = getNextGap(gap);

		// Initialize swapped as false so that we can
		// check if swap happened or not
		swapped = false;

		// Compare all elements with current gap
		for (int i = 0; i < n - gap; i++)
		{
			if (a[i] > a[i + gap])
			{
				std::swap(a[i], a[i + gap]);
				swapped = true;
			}
		}
	}
}

int binarySearch(int a[], int item, int low, int high)
{
	while (low <= high) {
		int mid = low + (high - low) / 2;
		if (item == a[mid])
			return mid + 1;
		else if (item > a[mid])
			low = mid + 1;
		else
			high = mid - 1;
	}

	return low;
}
void binInsertionSort(int a[], int n)
{
	int i, loc, j, k, selected;

	for (i = 1; i < n; ++i) {
		j = i - 1;
		selected = a[i];

		loc = binarySearch(a, selected, 0, j);

		while (j >= loc) {
			a[j + 1] = a[j];
			j--;
		}
		a[j + 1] = selected;
	}
}

int _findRandomPivot(int arr[], int start, int end)
{
	int n = end - start + 1;
	// Selecting the random pivot index
	int pivotInd = rand() % n;
	std::swap(arr[end], arr[start + pivotInd]);
	int pivot = arr[end];
	//initialising pivoting point to start index
	pivotInd = start;
	for (int i = start; i < end; i++) {

		// If an element is lesser than pivot, swap it.
		if (arr[i] <= pivot) {
			std::swap(arr[i], arr[pivotInd]);

			// Incrementing pivotIndex for further
			// swapping.
			pivotInd++;
		}
	}

	// Lastly swapping or the
	// correct position of pivot
	std::swap(arr[pivotInd], arr[end]);
	return pivotInd;
}
int _partition(int arr[], int low, int high)
{
	int pivot = arr[high]; //pivot is last elem
	//int pivot = findRandomPivot(arr, low, high);
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++)
	{
		if (arr[j] < pivot)
		{
			++i;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}
void _quickSort(int arr[], int low, int high)
{
	if (low < high)
	{
		int pi = _partition(arr, low, high);

		_quickSort(arr, low, pi - 1);
		_quickSort(arr, pi + 1, high);
	}
}

int partitionHoare(int arr[], int low, int high)
{
	int pivot = arr[low];
	int i = low - 1, j = high + 1;

	while (true) {

		// Find leftmost element greater than
		// or equal to pivot
		do {
			i++;
		} while (arr[i] < pivot);

		// Find rightmost element smaller than
		// or equal to pivot
		do {
			j--;
		} while (arr[j] > pivot);

		// If two pointers met
		if (i >= j)
			return j;

		std::swap(arr[i], arr[j]);
	}
}
int partition_rHoare(int arr[], int low, int high)
{
	int random = low + rand() % (high - low + 1);

	// Swap A[random] with A[high]
	std::swap(arr[random], arr[low]);

	return partitionHoare(arr, low, high);
}
int partition_mHoare(int arr[], int low, int high)
{
	//int elem = low + rand() % (high - low);
	int elem = (low + high) / 2;
	if ((arr[elem] > arr[high]) ^ (arr[elem] > arr[low]))
		std::swap(arr[low], arr[elem]);
	else if ((arr[high] < arr[elem]) ^ (arr[high] < arr[low]))
		std::swap(arr[low], arr[high]);
	/*
	if (arr[elem] < arr[low])
		std::swap(arr[low], arr[elem]);
	if (arr[high] < arr[low])
		std::swap(arr[low], arr[high]);
	if (arr[elem] < arr[high])
		std::swap(arr[elem], arr[high]);
	std::swap(arr[low], arr[high]);
	*/
	return partitionHoare(arr, low, high);
}
void quickSortHoare_r(int arr[], int low, int high)
{
	if (low < high) {
		// pi is partitioning index,
		// arr[p] is now at right place
		int pi = partition_rHoare(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quickSortHoare_r(arr, low, pi);
		quickSortHoare_r(arr, pi + 1, high);
	}
}
void quickSortHoare_m(int arr[], int low, int high)
{
	if (low < high) {
		// pi is partitioning index,
		// arr[p] is now at right place
		int pi = partition_mHoare(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quickSortHoare_m(arr, low, pi);
		quickSortHoare_m(arr, pi + 1, high);
	}
}
void quickSortHoare_f(int arr[], int low, int high)
{
	if (low < high) {
		// pi is partitioning index,
		// arr[p] is now at right place
		int pi = partitionHoare(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quickSortHoare_f(arr, low, pi);
		quickSortHoare_f(arr, pi + 1, high);
	}
}

int partitionLomuto(int arr[], int low, int high)
{
	int pivot = arr[high];    // pivot
	int i = (low - 1);  // Index of smaller element

	for (int j = low; j <= high - 1; j++)
	{
		// If current element is smaller than or
		// equal to pivot
		if (arr[j] <= pivot)
		{
			i++;    // increment index of smaller element
			std::swap(arr[i], arr[j]);
		}
	}
	std::swap(arr[i + 1], arr[high]);
	return (i + 1);
}
int partition_rLomuto(int arr[], int low, int high)
{
	int random = low + rand() % (high - low + 1);

	// Swap A[random] with A[high]
	std::swap(arr[random], arr[low]);

	return partitionLomuto(arr, low, high);
}
int partition_mLomuto(int arr[], int low, int high)
{
	//int elem = low + rand() % (high - low);
	int elem = (low + high) / 2;
	if (arr[elem] < arr[low])
		std::swap(arr[low], arr[elem]);
	if (arr[high] < arr[low])
		std::swap(arr[low], arr[high]);
	if (arr[elem] < arr[high])
		std::swap(arr[elem], arr[high]);

	return partitionLomuto(arr, low, high);
}
void quickSortLomuto_r(int arr[], int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
		   at right place */
		int pi = partition_rLomuto(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quickSortLomuto_r(arr, low, pi - 1);
		quickSortLomuto_r(arr, pi + 1, high);
	}
}
void quickSortLomuto_m(int arr[], int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
		   at right place */
		int pi = partition_mLomuto(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quickSortLomuto_m(arr, low, pi - 1);
		quickSortLomuto_m(arr, pi + 1, high);
	}
}
void quickSortLomuto_f(int arr[], int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
		   at right place */
		int pi = partitionLomuto(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		quickSortLomuto_f(arr, low, pi - 1);
		quickSortLomuto_f(arr, pi + 1, high);
	}
}

int larr[100000000], rarr[100000000];
void merge(int *array, int l, int m, int r) {
	int i, j, k, nl, nr;
	nl = m - l + 1; nr = r - m;
	//int *larr = new int[nl], *rarr = new int[nr];

	for (i = 0; i < nl; i++)
		larr[i] = array[l + i];
	for (j = 0; j < nr; j++)
		rarr[j] = array[m + 1 + j];

	i = 0; j = 0; k = l;
	while (i < nl && j < nr) {
		if (larr[i] <= rarr[j]) {
			array[k] = larr[i];
			i++;
		}
		else {
			array[k] = rarr[j];
			j++;
		}
		k++;
	}
	while (i < nl) {
		array[k] = larr[i];
		i++; k++;
	}
	while (j < nr) {
		array[k] = rarr[j];
		j++; k++;
	}

	//delete[] larr;
	//delete[] rarr;
}
void mergeSort(int *array, int l, int r) {
	if (l < r) {
		int m = l + (r - l) / 2;
		mergeSort(array, l, m);
		mergeSort(array, m + 1, r);
		merge(array, l, m, r);
	}
}

void mmerge(int *array, int l, int m, int r) {
	int i, j, k, nl, nr;
	nl = m - l + 1; nr = r - m;
	int *plarr = new int[nl], *prarr = new int[nr];

	for (i = 0; i < nl; i++)
		plarr[i] = array[l + i];
	for (j = 0; j < nr; j++)
		prarr[j] = array[m + 1 + j];

	i = 0; j = 0; k = l;
	while (i < nl && j < nr) {
		if (plarr[i] <= prarr[j]) {
			array[k] = plarr[i];
			i++;
		}
		else {
			array[k] = prarr[j];
			j++;
		}
		k++;
	}
	while (i < nl) {
		array[k] = plarr[i];
		i++; k++;
	}
	while (j < nr) {
		array[k] = prarr[j];
		j++; k++;
	}

	//delete[] larr;
	//delete[] rarr;
}
void mmergeSort(int *array, int l, int r) {
	if (l < r) {
		int m = l + (r - l) / 2;
		mmergeSort(array, l, m);
		mmergeSort(array, m + 1, r);
		mmerge(array, l, m, r);
	}
}

void it_merge(int *array, int l, int m, int r, int *plarr, int *prarr) {
	int i, j, k, nl, nr;
	nl = m - l + 1; nr = r - m;

	for (i = 0; i < nl; i++)
		plarr[i] = array[l + i];
	for (j = 0; j < nr; j++)
		prarr[j] = array[m + 1 + j];

	i = 0; j = 0; k = l;
	while (i < nl && j < nr) {
		if (plarr[i] <= prarr[j]) {
			array[k] = plarr[i];
			i++;
		}
		else {
			array[k] = prarr[j];
			j++;
		}
		k++;
	}
	while (i < nl) {
		array[k] = plarr[i];
		i++; k++;
	}
	while (j < nr) {
		array[k] = prarr[j];
		j++; k++;
	}

}
void iterativeMergeSort(int arr[], int n)
{
	int curr_size;
	int left_start;

	for (curr_size = 1; curr_size <= n - 1; curr_size = 2 * curr_size)
	{
		// Pick starting point of different subarrays of current size
		for (left_start = 0; left_start < n - 1; left_start += 2 * curr_size)
		{
			// Find ending point of left subarray. mid+1 is starting
			// point of right
			int mid = std::min(left_start + curr_size - 1, n - 1);

			int right_end = std::min(left_start + 2 * curr_size - 1, n - 1);

			// Merge Subarrays arr[left_start...mid] & arr[mid+1...right_end]
			merge(arr, left_start, mid, right_end);
		}
	}

	//delete[] plarr;
	//delete[] prarr;
}


void heapify(int arr[], int n, int i)
{
	int largest = i; // Initialize largest as root
	int l = 2 * i + 1; // left = 2*i + 1
	int r = 2 * i + 2; // right = 2*i + 2

	// If left child is larger than root
	if (l < n && arr[l] > arr[largest])
		largest = l;

	// If right child is larger than largest so far
	if (r < n && arr[r] > arr[largest])
		largest = r;

	// If largest is not root
	if (largest != i) {
		std::swap(arr[i], arr[largest]);

		// Recursively heapify the affected sub-tree
		heapify(arr, n, largest);
	}
}
void heapSort(int arr[], int n)
{
	// Build heap (rearrange array)
	for (int i = n / 2 - 1; i >= 0; i--)
		heapify(arr, n, i);

	// One by one extract an element from heap
	for (int i = n - 1; i > 0; i--) {
		// Move current root to end
		std::swap(arr[0], arr[i]);

		// call max heapify on the reduced heap
		heapify(arr, i, 0);
	}
}

int getMax(int arr[], int n)
{
	int mx = arr[0];
	for (int i = 1; i < n; i++)
		if (arr[i] > mx)
			mx = arr[i];
	return mx;
}
int output[10000000];
void countSort(int arr[], int n, int exp)
{
	int i, count[10] = { 0 };

	// Store count of occurrences in count[]
	for (i = 0; i < n; i++)
		count[(arr[i] / exp) % 10]++;

	// Change count[i] so that count[i] now contains actual
	//  position of this digit in output[]
	for (i = 1; i < 10; i++)
		count[i] += count[i - 1];

	// Build the output array
	for (i = n - 1; i >= 0; i--) {
		output[count[(arr[i] / exp) % 10] - 1] = arr[i];
		count[(arr[i] / exp) % 10]--;
	}

	// Copy the output array to arr[], so that arr[] now
	// contains sorted numbers according to current digit
	for (i = 0; i < n; i++)
		arr[i] = output[i];
}
void radixsort(int arr[], int n)
{
	// Find the maximum number to know number of digits
	int m = getMax(arr, n);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
		countSort(arr, n, exp);
}


void insertionSort(int arr[], int left, int right)
{
	for (int i = left + 1; i <= right; ++i)
	{
		int temp = arr[i];
		int j = i - 1;
		while (j >= left && arr[j] > temp)
		{
			arr[j + 1] = arr[j];
			--j;
		}
		arr[j + 1] = temp;
	}
}
void timSort(int arr[], int n, const int RUN)
{
	for (int i = 0; i < n; i += RUN)
		insertionSort(arr, i, std::min((i + RUN - 1), (n - 1)));

	for (int size = RUN; size < n; size = 2 * size)
	{
		for (int left = 0; left < n; left += 2 * size)
		{
			int mid = left + size - 1;
			int right = std::min((left + 2 * size - 1), (n - 1));

			if (mid < right)
				merge(arr, left, mid, right);
		}
	}
}

void radixtimSort(int arr[], int n, const int RUN)
{
	for (int i = 0; i < n; i += RUN)
		//insertionSort(arr, i, std::min((i + RUN - 1), (n - 1)));
		radixsort(arr + i, std::min(RUN, n - i));

	for (int size = RUN; size < n; size = 2 * size)
	{
		for (int left = 0; left < n; left += 2 * size)
		{
			int mid = left + size - 1;
			int right = std::min((left + 2 * size - 1), (n - 1));

			if (mid < right)
				merge(arr, left, mid, right);
		}
	}
}

struct Node
{
	int key;
	struct Node *left, *right;
};
struct Node *newNode(int item)
{
	struct Node *temp = new Node;
	temp->key = item;
	temp->left = temp->right = NULL;
	return temp;
}
void storeSorted(Node *root, int arr[], int &i)
{
	if (root != NULL)
	{
		storeSorted(root->left, arr, i);
		delete root->left;
		arr[i++] = root->key;
		storeSorted(root->right, arr, i);
		delete root->right;
	}
}
Node* insert(Node* node, int key)
{
	/* If the tree is empty, return a new Node */
	if (node == NULL) return newNode(key);

	/* Otherwise, recur down the tree */
	if (key < node->key)
		node->left = insert(node->left, key);
	else if (key > node->key)
		node->right = insert(node->right, key);

	/* return the (unchanged) Node pointer */
	return node;
}
void treeSort(int arr[], int n)
{
	struct Node *root = NULL;

	// Construct the BST
	root = insert(root, arr[0]);
	for (int i = 1; i < n; i++)
		root = insert(root, arr[i]);

	int i = 0;
	storeSorted(root, arr, i);
}

//new sort
void newSort1(int arr[], int low, int high, int MLEN)
{
	if (low < high - 1)
	{


		const int size = high - low + 1;

		if (size < MLEN)
		{

			insertionSort(arr, low, high);
		}

		else
		{
			int pi = partition_mHoare(arr, low, high);
			newSort1(arr, low, pi, MLEN);
			newSort1(arr, pi + 1, high, MLEN);
			//heapSort(arr + pi + 1, size);
		}
	}
}
void newSort2(int arr[], int low, int high, int depth)
{
	if (low < high - 1)
	{

		
		const int size = high - low + 1;

		if (size < 32)
		{

			insertionSort(arr, low, high);
		}
		else if (depth > 16)
			timSort(arr + low, size, 32);
		else 
		{
			int pi = partition_mHoare(arr, low, high);
			if (arr[pi] < 100)
				radixsort(arr + low, pi - low + 1);
			else
				newSort2(arr, low, pi, depth + 1);
				newSort2(arr, pi + 1, high, depth + 1);
		}
	}
}


int partition(int arr[], int low, int high)
{
	//median-of-three approach
	int elem = (low + high) / 2;
	if ((arr[elem] > arr[high]) ^ (arr[elem] > arr[low]))
		std::swap(arr[low], arr[elem]);
	else if ((arr[high] < arr[elem]) ^ (arr[high] < arr[low]))
		std::swap(arr[low], arr[high]);

	//partitioning
	int pivot = arr[low];
	int i = low, j = high;
	while (i < j) {

		//find leftmost element greater than or equal to the pivot
		while (arr[i] < pivot)
			++i;

		//find rightmost element smaller than or equal to the pivot
		while (arr[j] > pivot)
			--j;

		std::swap(arr[i], arr[j]); //swap the two elements
	}
	return j;
}
void sort(int arr[], const int low, const int high, const int depth, int MLEN)
{
	if (low < high - 1) //check if the length is greater than 2
	{
		const int size = high - low + 1; 

		if (size < MLEN)
			insertionSort(arr, low, high); //do Insertion Sort on small arrays

		else if (depth == 0);
			//heapSort(arr, size); //do Heap Sort when the maximum depth is exceeded

		else //do Quicksort otherwise, keeping count of the depth
		{
			const int p = partition(arr, low, high);
			sort(arr, low, p, depth - 1, MLEN);
			sort(arr, p + 1, high, depth - 1, MLEN);
		}
	}
}

void sort(int arr[], int n)
{
	//sort sub-arrays with Insertion Sort
	for (int i = 0; i < n; i += 32)
		insertionSort(arr, i, std::min((i + 32 - 1), (n - 1)));

	//merge the sub-arrays
	for (int size = 32; size < n; size = size << 1)
	{
		for (int low = 0; low < n; low += size << 1)
		{
			int mid = low + size - 1, high = std::min((low + 2 * size - 1), (n - 1));

			if (mid < high)
				merge(arr, low, mid, high);
		}
	}
}


int main()
{
	srand(time(NULL));

	sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "CPlot");
	window.setFramerateLimit(30);

	sf::RectangleShape background(sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHT));
	background.setFillColor(sf::Color::White);


	Graph g(sf::Vector2f(100.0, 70.0), sf::Vector2f(1200.0, 900.0), 100000000, 7000000, 0, 0);
	g.addPlot("Iterative merge sort (static)", sf::Color::Red);
	g.addPlot("Hoare Quicksort (median-of-3)", sf::Color(225, 128, 0, 255));
	g.addPlot("Introsort (32 max size)", sf::Color(220, 0, 220, 255));
	g.addPlot("Timsort (32 element block)", sf::Color(0, 180, 0, 255));
	g.addPlot("", sf::Color(80, 230, 80, 255));
	g.addPlot("", sf::Color::Cyan);
	g.addPlot("", sf::Color::Yellow);



	const int size = 100000000, tests = 2;

	int arr[tests][size + 1], copyArr[tests][size + 1];
	for (int i = 0; i < size + 1; i += 10000000)
	{
		for (int j = 0; j < tests; j++)
		{
			genRandArray(arr[j], i);
			copyArray(copyArr[j], arr[j], i);
		}
		std::cout << "-----------\n";

		sf::Clock clock;
		sf::Time start;

		std::cout << " - Introsort (8): ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
		{
			//printArray(arr[j], i);
			//timSort(arr[j], i, 32);
			iterativeMergeSort(arr[j], i);
			//timSort(arr[j], i, 32);
			//std::sort(arr[j], arr[j] + i);
			//newSort1(arr[j], 0, i - 1, 8);
			//printArray(arr[j], i);
		}
		g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 0);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);

		std::cout << " - Introsort (16): ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
		{
			//printArray(arr[j], i);
			//newSort1(arr[j], 0, i - 1, 16);
			//timSort(arr[j], i, 32);
			quickSortHoare_m(arr[j], 0, i - 1);
			//printArray(arr[j], i);
		}
		g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 1);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Introsort (32): ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
		{
			//printArray(arr[j], i);
			//radixsort(arr[j], i);
			//heapSort(arr[j], i);
			//std::sort(arr[j], arr[j] + i);
			newSort1(arr[j], 0, i - 1, 32);
			//insertionSort(arr[j], i);
			//printArray(arr[j], i);
		}
		g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 2);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Introsort (64): ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
		{
			//mergeSort(arr[j], 0, i - 1);
			//quickSortHoare_m(arr[j], 0, i - 1);
			//newSort1(arr[j], 0, i - 1, 64);
			//shellSort(arr[j], i);
			timSort(arr[j], i, 32);
			//printArray(arr[j], i);

		}
		g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 3);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Introsort (128): ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
		{
			//printArray(arr[j], i);
			//mergeSort(arr[j], 0, i - 1);
			//quickSortLomuto_m(arr[j], 0, i - 1);
			//iterativeMergeSort(arr[j], i);
			//newSort2(arr[j], 0, i - 1, 0);
			//newSort1(arr[j], 0, i - 1, 128);
			//printArray(arr[j], i);
		}
		//g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 4);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Heap Sort: ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
		{
			//quickSortHoare_m(arr[j], 0, i - 1);
		}
		//g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 5);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Binary Insertion Sort: ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j);
		//binInsertionSort(arr[j], i);
	//g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 6);
		std::cout << "DONE\n";

		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Shell sort: ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j);
		//insertionSort(arr[j], i);
	//g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 6);
		std::cout << "DONE\n";

		/*
		for (int j = 0; j < tests; ++j)
			copyArray(arr[j], copyArr[j], i);
		std::cout << " - Gnome sort: ";
		start = clock.getElapsedTime();
		for (int j = 0; j < tests; ++j)
			gnomeSort(arr[j], i);
		g.addPointToPlot(sf::Vector2f(i, (clock.getElapsedTime().asMicroseconds() - start.asMicroseconds()) / tests), 8);
		std::cout << "DONE\n";
		*/
	}



	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear();
		window.draw(background);
		g.draw(window);

		window.display();
	}

	return 0;
}

/*
Ideas:
	- use a simple O(n) alg to sort lenght 3 subarrays followed by insertion/merge
	- is almost sorted find the usorted subarray
	- compute a sort_ratio in O(n)
	- if reversed, reverse


Ignored:
	- count (too specific)
	- bucket (radix for floats)
	- cycle (if swap is costly, not our case)
	- strand (not in place)
	- bitonic (no reason in a single thread)

For O(n*n) compare:
  - selection
  - bubble
  - comb (improv bubble)
  - cocktail (improv bubble 2)
  - odd-even (improv bubble 3)
  - insertion
  - shell (improv insert)
  - bin insertion (improv insert)
  - gnome (bad but...)

  For O(n*logn) compare:
  - quick (fixed pivot)
  - quick (rand pivot)
  - quick (median-of-three pivot + median-of-three rand)
  - quick (median pivot)
  - quick (hoare/lomuto pivot)
  - merge
  - iterative merge
  - heap
  - tree

  For hybrid:
  - tim (merge + insert)

  -radix sort
 */