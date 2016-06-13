#ifndef COMPARABLE_H__
#define COMPARABLE_H__

template<typename T>
class Comparable { public:
	Bit operator>(const T&rhs) const {
		return static_cast<const T*>(this)->greater(rhs);
	}
	Bit operator<(const T& rhs) const {
		return rhs > (*static_cast<const T*>(this));
	}

	Bit operator<=(const T& rhs) const {
		return !(*static_cast<const T*>(this) > rhs); 
	}

	Bit operator>=(const T& rhs) const {
		return !(rhs > *static_cast<const T*>(this));
	}
	
	Bit operator==(const T& rhs) const {
		return static_cast<const T*>(this)->equal(rhs);
	}
	Bit operator!=(const T& rhs) const {
		return !(*static_cast<const T*>(this) == rhs);
	}
};
#endif