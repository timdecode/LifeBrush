//
//  dbscan.hpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2016-01-01.
//  Copyright © 2016 Timothy Davison, Inc. All rights reserved.
//

#pragma once

#include <vector>
#include <unordered_set>

#include "PointCloud.hpp"

namespace rg
{
	template<typename Scalar, typename Point>
	struct DBScan
	{
		struct PointAndIndex
		{
			PointAndIndex(Point p, size_t i) : point(p), index(i) {}

			Point point;
			size_t index;
		};

		static void dbscan(Point_kdTree<Scalar, Point>& index,
			Scalar radius,
			int minPoints,
			std::vector<std::vector<PointAndIndex>>& clusters_out,
			std::vector<PointAndIndex>& noise_out)
		{
			clusters_out.clear();
			noise_out.clear();

			auto& points = index.cloud.pts;
			auto& kdTree = index.kdTree;

			if (points.size() == 0)
				return;

			float radiusSqrd = radius * radius;

			nanoflann::SearchParams searchParams;
			searchParams.sorted = false;

			std::unordered_set<size_t> unvisited;
			for (int i = 0; i < points.size(); ++i)
				unvisited.insert(i);

			std::unordered_set<size_t> clustered;

			while (unvisited.size())
			{
				size_t i = *unvisited.begin();
				unvisited.erase(i);

				Point& p = points[i];

				std::vector<std::pair<size_t, float> > neighbours;
				kdTree.radiusSearch(&p[0], radiusSqrd, neighbours, searchParams);

				// -1 because we include the query point
				if (neighbours.size() - 1 < minPoints)
					noise_out.emplace_back(p, i);
				else
				{
					clusters_out.emplace_back();
					auto& cluster = clusters_out.back();

					// expand cluster
					cluster.emplace_back(p, i);
					clustered.insert(i);

					for (auto& pair : neighbours)
					{
						size_t i_ = pair.first;
						Point& p_ = points[i_];

						if (unvisited.find(i_) != unvisited.end())
						{
							unvisited.erase(i_);

							std::vector<std::pair<size_t, float> > neighbours_;
							kdTree.radiusSearch(&p_[0], radiusSqrd, neighbours, searchParams);

							// -1 because we include the query point
							if (neighbours_.size() - 1 >= minPoints)
								neighbours.insert(neighbours.end(), neighbours_.begin(), neighbours_.end());

							if (clustered.find(i_) == clustered.end())
								cluster.emplace_back(p_, i_);
						}
					}
				}
			}
		}
	};


};