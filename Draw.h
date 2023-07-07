#pragma once
#include"Cal.h"
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidht = 2.0f;
	const int kSoubdibision = 10;
	const float kGridEvery = (kGridHalfWidht * 2.0f) / float(kSoubdibision);

	Vector3 worldVertev[2];

	Vector3 screenVertex[2];
	Vector3 ndcVertex;
	for (int xIndex = 0; xIndex <= kSoubdibision; xIndex++)
	{
		worldVertev[0] = { xIndex * kGridEvery - kGridHalfWidht, 0.0f , -kGridHalfWidht };
		worldVertev[1] = { xIndex * kGridEvery - kGridHalfWidht, 0.0f, kGridHalfWidht };

		for (int i = 0; i < 2; i++)
		{
			ndcVertex = Transform(worldVertev[i], viewProjectionMatrix);
			screenVertex[i] = Transform(ndcVertex, viewportMatrix);
		}
		Novice::DrawLine(int(screenVertex[0].x), int(screenVertex[0].y), int(screenVertex[1].x), int(screenVertex[1].y), 0xAAAAAAFF);
	}

	for (int zIndex = 0; zIndex <= kSoubdibision; zIndex++)
	{
		worldVertev[0] = { -kGridHalfWidht, 0 , zIndex * kGridEvery - kGridHalfWidht };
		worldVertev[1] = { kGridHalfWidht, 0 , zIndex * kGridEvery - kGridHalfWidht };

		for (int i = 0; i < 2; i++)
		{
			ndcVertex = Transform(worldVertev[i], viewProjectionMatrix);
			screenVertex[i] = Transform(ndcVertex, viewportMatrix);
		}
		Novice::DrawLine(int(screenVertex[0].x), int(screenVertex[0].y), int(screenVertex[1].x), int(screenVertex[1].y), 0xAAAAAAFF);
	}
}

void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color)
{
	const uint32_t kSubdivision = 30;
	const float kLonEvery = 2.0f * float(std::numbers::pi) / float(kSubdivision);//経度分割1つ分の角度
	const float kLatEvery = float(std::numbers::pi) / float(kSubdivision);//緯度分割1つ分の角度
	//緯度の方向に分割 
	for (uint32_t latIndex = 0; latIndex < kSubdivision; latIndex++) {
		float lat = float(-std::numbers::pi) / 2.0f + kLatEvery * latIndex;//現在の緯度
		//経度の方向に分割
		for (uint32_t lonIndex = 0; lonIndex < kSubdivision; lonIndex++) {
			float lon = lonIndex * kLonEvery;//現在の経度
			//world座標系でのa,b,cを求める
			Vector3 a, b, c;
			a = { sphere.radius * std::cosf(lat) * std::cosf(lon), sphere.radius * std::sinf(lat), sphere.radius * std::cosf(lat) * std::sinf(lon) };
			a = Add(a, sphere.center);
			b = { sphere.radius * std::cosf(lat + kLatEvery) * std::cosf(lon), sphere.radius * std::sinf(lat + kLatEvery), sphere.radius * std::cosf(lat + kLatEvery) * std::sinf(lon) };
			b = Add(b, sphere.center);
			c = { sphere.radius * std::cosf(lat) * std::cosf(lon + kLonEvery), sphere.radius * std::sinf(lat), sphere.radius * std::cosf(lat) * std::sinf(lon + kLonEvery) };
			c = Add(c, sphere.center);

			//a,b,cをスクリーン座標へ
			a = Transform(a, viewProjectionMatrix);
			a = Transform(a, viewportMatrix);
			b = Transform(b, viewProjectionMatrix);
			b = Transform(b, viewportMatrix);
			c = Transform(c, viewProjectionMatrix);
			c = Transform(c, viewportMatrix);


			//線を引く
			Novice::DrawLine(
				int(a.x), int(a.y),
				int(b.x), int(b.y),
				color
			);

			Novice::DrawLine(
				int(a.x), int(a.y),
				int(c.x), int(c.y),
				color
			);
		}
	}
};

void DrawPlane(const Plane& plane, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix, uint32_t color) {

	Vector3 center = F2VMultiply(plane.distance, plane.normal);
	Vector3 perpendiculars[4];
	perpendiculars[0] = Normalize(Prependicular(plane.normal));
	perpendiculars[1] = { -perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z };
	perpendiculars[2] = Cross(plane.normal, perpendiculars[0]);
	perpendiculars[3] = { -perpendiculars[2].x, -perpendiculars[2].y, -perpendiculars[2].z };

	Vector3 points[4];
	for (uint32_t index = 0; index < 4; ++index) {
		Vector3 extend = F2VMultiply(2.0f, perpendiculars[index]);
		Vector3 point = Add(center, extend);
		points[index] = Transform(Transform(point, viewProjectionMatrix), viewportMatrix);
	}

	Novice::DrawLine(int(points[0].x), int(points[0].y), int(points[3].x), int(points[3].y), color);
	Novice::DrawLine(int(points[3].x), int(points[3].y), int(points[1].x), int(points[1].y), color);
	Novice::DrawLine(int(points[1].x), int(points[1].y), int(points[2].x), int(points[2].y), color);
	Novice::DrawLine(int(points[2].x), int(points[2].y), int(points[0].x), int(points[0].y), color);

};