#ifndef THREAD_HPP_
#define THREAD_HPP_
/*
 * This implementation is copied from
 * http://stackoverflow.com/questions/12215395/thread-pool-using-boost-asio
 *
 */
#include <boost/asio.hpp>
#include <boost/thread.hpp>

class thread_pool
{
private:
  boost::asio::io_service io_service_;
  boost::asio::io_service::work work_;
  boost::thread_group threads_;
  std::size_t available_;
  boost::mutex mutex_;
public:

  /// @brief Constructor.
  thread_pool( std::size_t pool_size )
    : work_( io_service_ ),
      available_( pool_size )
  {
    for ( std::size_t i = 0; i < pool_size; ++i )
    {
      threads_.create_thread( boost::bind( &boost::asio::io_service::run,
                                           &io_service_ ) );
    }
  }

  /// @brief Destructor.
  ~thread_pool()
  {
    // Force all threads to return from io_service::run().
    io_service_.stop();

    // Suppress all exceptions.
    try
    {
      threads_.join_all();
    }
    catch ( ... ) {}
  }

  /// @brief Adds a task to the thread pool if a thread is currently available.
  template < typename Task >
  void run_task( Task task )
  {
    boost::unique_lock< boost::mutex > lock( mutex_ );

    // If no threads are available, then return.
    if ( 0 == available_ ) return;

    // Decrement count, indicating thread is no longer available.
    --available_;

    // Post a wrapped task into the queue.
    io_service_.post( boost::bind( &thread_pool::wrap_task, this,
                                   boost::function< void() >( task ) ) );
  }

private:
  /// @brief Wrap a task so that the available count can be increased once
  ///        the user provided task has completed.
  void wrap_task( boost::function< void() > task )
  {
    // Run the user supplied task.
    try
    {
      task();
    }
    // Suppress all exceptions.
    catch ( ... ) {}

    // Task has finished, so increment count of available threads.
    boost::unique_lock< boost::mutex > lock( mutex_ );
    ++available_;
  }
};

#endif /* THREAD_HPP_ */
