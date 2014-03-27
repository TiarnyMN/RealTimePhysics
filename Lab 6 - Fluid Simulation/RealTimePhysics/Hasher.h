class Hasher
{
public:
	int operator()(glm::vec3 val)
	{
		int x = val.x * 73856093;
		int y = val.y * 19349663;
		int z = val.z * 83492791;
		
		return (x ^ y ^ z);
	}

	bool operator()(glm::vec3 first, glm::vec3 second)
	{
		return first == second;
	}
};