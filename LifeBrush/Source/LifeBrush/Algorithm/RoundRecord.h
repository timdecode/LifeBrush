// Copyright (c) 2017 Timothy Davison. All rights reserved.

#pragma once

#include <vector>
#include <memory>
#include <ctime>

class RoundRecord
{
public:
	std::clock_t start = 0.0f;
	std::clock_t end = 0.0f;

	auto duration() -> double
	{
		return (end - start) / (double)CLOCKS_PER_SEC;
	}


	float totalEnergy_kCoherence = 0.0f;
	float totalEnergy_bruteForce = 0.0f;

	float averageEnergy_kCoherence = 0.0f;
	float averageEnergy_bruteForce = 0.0f;

	uint32_t elementsTotal = 0;

	uint32_t elementsGenerated = 0;
	uint32_t elementsModified = 0;
	uint32_t elementsRemoved = 0;
};

class RecordSummary
{
public:
	virtual ~RecordSummary() {}

	virtual auto statsString() -> FString
	{
		FString result = "---";

		size_t ri = 0;
		for(auto& record : records)
		{
			result += "\n";

			result += FString::Printf( TEXT( "Round %i" ), ri );
			result += FString::Printf( TEXT( "  duration: %f\n" ), record->duration() );

			result += FString::Printf( TEXT( "  totalEnergy_kCoherence: %f\n" ), record->totalEnergy_kCoherence );
			result += FString::Printf( TEXT( "  averageEnergy_kCoherence: %f\n" ), record->averageEnergy_kCoherence );

			result += FString::Printf( TEXT( "  totalEnergy_bruteForce: %f\n" ), record->totalEnergy_bruteForce );
			result += FString::Printf( TEXT( "  averageEnergy_bruteForce: %f\n" ), record->averageEnergy_bruteForce );

			result += FString::Printf( TEXT( "  elementsTotal: %i\n" ), record->elementsTotal );

			result += FString::Printf( TEXT( "  elementsGenerated: %i\n" ), record->elementsGenerated );
			result += FString::Printf( TEXT( "  elementsModified: %i\n" ), record->elementsModified );
			result += FString::Printf( TEXT( "  elementsRemoved: %i\n" ), record->elementsRemoved );

			ri++;
		}

		result += FString::Printf( TEXT( "Total Duration: %f\n" ), duration() );

		return result;
	}

	virtual auto separatedValueStatsString(FString separator) -> FString
	{
		// concatentation, but its easier to read in source
		FString result = "round" + separator;
		result += "duration" + separator;

		result += "totalEnergy_kCoherence" + separator;
		result += "averageEnergy_kCoherence" + separator;

		result += "totalEnergy_bruteForce" + separator;
		result += "averageEnergy_bruteForce" + separator;

		result += "elementsTotal" + separator;
		result += "elementsGenerated" + separator;
		result += "elementsModified" + separator;
		result += "elementsRemoved" + separator;

		size_t ri = 0;
		for(auto& record : records)
		{
			result += "\n";

			result += FString::Printf( TEXT( "%i%s" ), ri, *separator );
			result += FString::Printf( TEXT( "%f%s" ), record->duration(), *separator );

			result += FString::Printf( TEXT( "%f%s" ), record->totalEnergy_kCoherence, *separator );
			result += FString::Printf( TEXT( "%f%s" ), record->averageEnergy_kCoherence, *separator );

			result += FString::Printf( TEXT( "%f%s" ), record->totalEnergy_bruteForce, *separator );
			result += FString::Printf( TEXT( "%f%s" ), record->averageEnergy_bruteForce, *separator );

			result += FString::Printf( TEXT( "%i%s" ), record->elementsTotal, *separator );

			result += FString::Printf( TEXT( "%i%s" ), record->elementsGenerated, *separator );
			result += FString::Printf( TEXT( "%i%s" ), record->elementsModified, *separator );
			result += FString::Printf( TEXT( "%i" ), record->elementsRemoved );

			ri++;
		}

		result += FString::Printf( TEXT( "\nTotal Duration: %f\n" ), duration() );

		return result;
	}

	virtual auto emplace() -> RoundRecord&
	{
		records.emplace_back( new RoundRecord() );

		return *(records.back().get());
	}

	virtual auto back() -> RoundRecord&
	{
		return *(records.back().get());
	}

	virtual void clear()
	{
		records.clear();
	}

	auto duration() -> double
	{
		if(records.size() == 0)
			return 0.0;

		double d = 0.0;

		for(auto& record : records)
		{
			d += record->duration();
		}

		return d;
	}

private:
	std::vector< std::unique_ptr<RoundRecord> > records;
};

class PatchRecordSummary : public RecordSummary
{
public:
	virtual ~PatchRecordSummary() {}
};