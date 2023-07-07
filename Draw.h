#pragma once
#include"Cal.h"
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;
	const uint32_t kSubdivision = 10;
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);

	Vector3 localBorderVer[2]{};
	Vector3 localStripeVer[2]{};

	/*for (uint32_t x = 0; x <= kSubdivision; x++) {
		localBorderVer[x][0] = { -kGridHalfWidth, 0.0f, kGridEvery * (float(x) - 5) };
		localBorderVer[x][1] = { kGridHalfWidth, 0.0f, kGridEvery * (float(x) - 5) };

		localStripeVer[x][0] = { kGridEvery * (float(x) - 5) , 0.0f, -kGridHalfWidth };
		localStripeVer[x][1] = { kGridEvery * (float(x) - 5) , 0.0f, kGridHalfWidth };

	}*/

	Vector3 screenBorderVer[2]{};
	Vector3 screenStripeVer[2]{};


	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {


		localBorderVer[0] = { -kGridHalfWidth, 0.0f, kGridEvery * (float(xIndex) - 5) };
		localBorderVer[1] = { kGridHalfWidth, 0.0f, kGridEvery * (float(xIndex) - 5) };

		localStripeVer[0] = { kGridEvery * (int(xIndex) - 5) , 0.0f, -kGridHalfWidth };
		localStripeVer[1] = { kGridEvery * (int(xIndex) - 5) , 0.0f, kGridHalfWidth };


		Vector3 ndcBorderStart = TransformCoord(localBorderVer[0], viewProjectionMatrix);
		Vector3 ndcBorderEnd = TransformCoord(localBorderVer[1], viewProjectionMatrix);

		Vector3 ndcStripeStart = TransformCoord(localStripeVer[0], viewProjectionMatrix);
		Vector3 ndcStripeEnd = TransformCoord(localStripeVer[1], viewProjectionMatrix);




		screenBorderVer[0] = TransformCoord(ndcBorderStart, viewportMatrix);
		screenBorderVer[1] = TransformCoord(ndcBorderEnd, viewportMatrix);

		screenStripeVer[0] = TransformCoord(ndcStripeStart, viewportMatrix);
		screenStripeVer[1] = TransformCoord(ndcStripeEnd, viewportMatrix);

		Novice::DrawLine(
			int(screenBorderVer[0].x), int(screenBorderVer[0].y),
			int(screenBorderVer[1].x), int(screenBorderVer[1].y),
			0xAAAAAAFF);

		Novice::DrawLine(
			int(screenStripeVer[0].x), int(screenStripeVer[0].y),
			int(screenStripeVer[1].x), int(screenStripeVer[1].y),
			0xAAAAAAFF);

		if (localBorderVer[0].z == 0) {
			Novice::DrawLine(
				int(screenStripeVer[0].x), int(screenStripeVer[0].y),
				int(screenStripeVer[1].x), int(screenStripeVer[1].y),
				0x000000FF);

			Novice::DrawLine(
				int(screenBorderVer[0].x), int(screenBorderVer[0].y),
				int(screenBorderVer[1].x), int(screenBorderVer[1].y),
				0x000000FF);
		}

	}

};

void DrawShere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	const uint32_t kSubdivision = 10;
	const float kLonEvery = (2 * 3.14f) / kSubdivision;
	const float kLatEvery = 3.14f / kSubdivision;
	//
	for (uint32_t latIndex = 0; latIndex < kSubdivision; ++latIndex) {
		float lat = -3.14f / 2.0f + kLatEvery * latIndex;
		//
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery;
			//
			Vector3 a, b, c;
			a.x = sphere.center.x + sphere.radius * std::cosf(lat) * std::cosf(lon);
			a.y = sphere.center.y + sphere.radius * std::sinf(lat);
			a.z = sphere.center.z + sphere.radius * std::cosf(lat) * std::sinf(lon);

			b.x = sphere.center.x + sphere.radius * std::cosf(lat + (3.14f / kSubdivision)) * std::cosf(lon);
			b.y = sphere.center.y + sphere.radius * std::sinf(lat + (3.14f / kSubdivision));
			b.z = sphere.center.z + sphere.radius * std::cosf(lat + (3.14f / kSubdivision)) * std::sinf(lon);

			c.x = sphere.center.x + sphere.radius * std::cosf(lat) * std::cosf(lon + ((3.14f * 2) / kSubdivision));
			c.y = sphere.center.y + sphere.radius * std::sinf(lat);
			c.z = sphere.center.z + sphere.radius * std::cosf(lat) * std::sinf(lon + ((3.14f * 2) / kSubdivision));

			Vector3 ndcA = TransformCoord(a, viewProjectionMatrix);
			Vector3 ndcB = TransformCoord(b, viewProjectionMatrix);
			Vector3 ndcC = TransformCoord(c, viewProjectionMatrix);

			Vector3 screenA = TransformCoord(ndcA, viewportMatrix);
			Vector3 screenB = TransformCoord(ndcB, viewportMatrix);
			Vector3 screenC = TransformCoord(ndcC, viewportMatrix);

			//
			//
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenB.x), int(screenB.y), color);
			Novice::DrawLine(int(screenA.x), int(screenA.y), int(screenC.x), int(screenC.y), color);

		}
	}
};
void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {
	Vector3 center = Multiply(plane.distance, plane.normal);
	Vector3 perpendiculars[4];
	perpendiculars[0] = Normalize(Perpendicular(plane.normal));
	perpendiculars[1] = { -perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z };
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
	perpendiculars[3] = { -perpendiculars[2].x, -perpendiculars[2].y, -perpendiculars[2].z };

	Vector3 points[4];
	for (uint32_t index = 0; index < 4; ++index) {
		Vector3 extend = Multiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = TransformCoord(TransformCoord(point, viewProjectionMatrix), viewportMatrix);
	}

	Novice::DrawLine(int(points[0].x), int(points[0].y), int(points[3].x), int(points[3].y), color);
	Novice::DrawLine(int(points[3].x), int(points[3].y), int(points[1].x), int(points[1].y), color);
	Novice::DrawLine(int(points[1].x), int(points[1].y), int(points[2].x), int(points[2].y), color);
	Novice::DrawLine(int(points[2].x), int(points[2].y), int(points[0].x), int(points[0].y), color);

};