// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

// critically damped positional spring
// see: http://mathproofs.blogspot.com/2013/07/critically-damped-spring-smoothing.html
struct PhysicsUtilities
{
	static FVector criticallyDampedSpringVelocity(
		const FVector v, 
		const FVector p, 
		const FVector target, 
		float w, 
		float dt)
	{
		return (v - std::pow(w, 2.0f) * dt * (p - target)) / std::pow(1 + w * dt, 2.0f);
	}

	static FVector alignmentTorque(const FQuat &dq, float& torqueMagnitude_out)
	{
		FVector torqueAxis;

		dq.ToAxisAndAngle(torqueAxis, torqueMagnitude_out);

		torqueMagnitude_out = std::sin(torqueMagnitude_out);

		return torqueAxis * torqueMagnitude_out;
	}
};