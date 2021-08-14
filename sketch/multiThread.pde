class MultiThreading extends Thread {
  public void run() {
    try {
      println("Thread " + Thread.currentThread().getId()+ " is running");
      //The stuff for the thread to process
    }
    catch (Exception e) {
      println("Exception is caught, stopping thread");
    }
  }
}

public class Multithread {
  ArrayList<MultiThreading> threads = new ArrayList<MultiThreading>();
  
  public void main(int nThreads) {
    for (int i = 0; i < nThreads; i++) {
      MultiThreading thread = new MultiThreading();
      threads.add(thread);
      thread.start();
    }
  }
  
  public void stopThread(int threadN) {
    threads.get(threadN).interrupt();
  }
}
