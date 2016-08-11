#define BOOST_TEST_MODULE dng::multithread

#include <dng/multithread.h>
#include <atomic>
#include <iostream>

BOOST_AUTO_TEST_CASE(test_basicpool) {
	using namespace std;
	using namespace dng::multithread;

	atomic<int> result{0};

	auto f = [&result](int x) {
		result += x;
	};

	{
		BasicPool<int> pool(f,4);
		for(int i=0;i<=100;++i) {
			pool.Enqueue(i);
		}
	}
	BOOST_CHECK(result == (100*101)/2);
}