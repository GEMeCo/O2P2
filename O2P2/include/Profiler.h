// ================================================================================================
// Basic instrumentation profiler by Yan Chernikov and David Churchill
//
// Singleton usage:
//
// // Begin session 
// Instrumentor::Get().BeginSession("Session Name");
// {
//	 // Place code like this in scopes you'd like to include in profiling:
//	 InstrumentationTimer timer("Profiled Scope Name");
//	 
//	 // or, if you want to use macros:
//	 PROFILE_FUNCTION()
//
//	 // Your code
// }
// // End Session
// Instrumentor::Get().EndSession();
//
// How to profile:
// This profiler will create a file named "results.json" only in debug mode.
// Just open the google chrome in address chrome://tracing and drag/paste the file.
// Pretty neat!
// 
// ================================================================================================
#pragma once

#include <string>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <thread>
#include <mutex>

#ifdef _DEBUG
#define PROFILE_SCOPE(name) InstrumentationTimer timer##__LINE__(name)
#define PROFILE_FUNCTION()  PROFILE_SCOPE(__FUNCTION__)
#else
#define PROFILE_SCOPE(name)
#define PROFILE_FUNCTION()
#endif

struct ProfileResult
{
	const std::string name;
	long long start, end;
	uint32_t threadID;
};

class Instrumentor
{
	std::string		m_sessionName = "None";
	std::ofstream	m_outputStream;
	int				m_profileCount = 0;
	std::mutex		m_lock;
	bool			m_activeSession = false;

public:

	static Instrumentor& instance() {
		static Instrumentor instance;
		return instance;
	};

	~Instrumentor() {
		endSession();
	};

	static void beginSession(const std::string& name, const std::string& filepath = "results.json") {
		instance().beginSessionImpl(name, filepath);
	};

	static void endSession() {
		instance().endSessionImpl();
	};

	static void writeProfile(const ProfileResult& result) {
		instance().writeProfileImpl(result);
	};

private:
	Instrumentor() { }

	void beginSessionImpl(const std::string& name, const std::string& filepath = "results.json") {
		if (m_activeSession) { endSession(); }
		m_activeSession = true;
		m_outputStream.open(filepath);
		writeHeader();
		m_sessionName = name;
	};

	void endSessionImpl() {
		if (!m_activeSession) { return; }
		m_activeSession = false;
		writeFooter();
		m_outputStream.close();
		m_profileCount = 0;
	};

	void writeProfileImpl(const ProfileResult& result) {
		std::lock_guard<std::mutex> lock(m_lock);

		if (m_profileCount++ > 0) { m_outputStream << ","; }

		std::string name = result.name;
		std::replace(name.begin(), name.end(), '"', '\'');

		m_outputStream << "{";
		m_outputStream << "\"cat\":\"function\",";
		m_outputStream << "\"dur\":" << (result.end - result.start) << ',';
		m_outputStream << "\"name\":\"" << name << "\",";
		m_outputStream << "\"ph\":\"X\",";
		m_outputStream << "\"pid\":0,";
		m_outputStream << "\"tid\":" << result.threadID << ",";
		m_outputStream << "\"ts\":" << result.start;
		m_outputStream << "}";
	};

	void writeHeader() {
		m_outputStream << "{\"otherData\": {},\"traceEvents\":[";
	};

	void writeFooter() {
		m_outputStream << "]}";
	};
};

class InstrumentationTimer
{
	ProfileResult m_result;

	std::chrono::time_point<std::chrono::high_resolution_clock> m_startTimepoint;
	bool m_stopped;

public:

	InstrumentationTimer(const std::string& name)
		: m_result({ name, 0, 0, 0 })
		, m_stopped(false)
	{
		m_startTimepoint = std::chrono::high_resolution_clock::now();
	}

	~InstrumentationTimer()
	{
		if (!m_stopped) { stop(); }
	}

	void stop()
	{
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		m_result.start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTimepoint).time_since_epoch().count();
		m_result.end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
		m_result.threadID = std::hash<std::thread::id>{}(std::this_thread::get_id());
		Instrumentor::writeProfile(m_result);

		m_stopped = true;
	}
};