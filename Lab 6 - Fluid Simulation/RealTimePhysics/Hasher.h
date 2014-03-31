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

struct Capsule
{
	glm::vec3 startPoint;
	glm::vec3 endPoint;

	float length;
	glm::vec3 centre;
	float radius;

	glm::vec3 q;

	glm::vec3 rotation;

	void CalculateCapsulePosition(glm::vec3 rot)
	{
		rotation += rot;
		startPoint = centre - glm::vec3(glm::sin(rotation.x) * length, glm::cos(rotation.x) * length, 0.0f);
		endPoint = centre + glm::vec3(glm::sin(rotation.x) * length, glm::cos(rotation.x) * length, 0.0f);
	}
};