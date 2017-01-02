class slow_index;
{
	private:
	size_t _1, _2, _3, _4;
	public:
	slow_index(const size_t &_N)
	{	
		_1=3*_N*_N*_N*2;
		_2=3*_N*_N*2;
		_3=3*_N*2;
		_4=3;
	}

	const size_t 	
	operator (const uint32_t &d, const uint32_t &X, const uint32_t &Y, const uint32_t &i, const uint32_t &j)
	return X*_1+Y*_2+d*_3+i*_4+j;
}

fast_bit_sum8
{
}

fast_bit_sum16
{
}

fast_bit_sum32
{
}

fast_bit_sum64
{
}

fast_bit_sum128
{
}

fast_bit_sumN
{
}

fast_index(const uint &d, const uint32_t &i, const uint32_t &j)
{
	return _index[d][i][j];
}
