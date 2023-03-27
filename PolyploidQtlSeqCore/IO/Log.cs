using Cysharp.Diagnostics;

namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// ログ
    /// </summary>
    internal static class Log
    {
        private static readonly List<string> _logs = new();
        private static readonly string _sepalator = new('#', 80);
        private static readonly object _syncObj = new();

        /// <summary>
        /// メッセージをログに追加する。
        /// </summary>
        /// <param name="message">メッセージ</param>
        public static void Add(string message)
        {
            lock (_syncObj) _logs.Add(message);
        }

        /// <summary>
        /// プロセスエラーログを追加する。
        /// </summary>
        /// <param name="ex">例外</param>
        public static void Add(ProcessErrorException ex)
        {
            lock (_syncObj)
            {
                _logs.AddRange(ex.ErrorOutput);
                _logs.Add(ex.ToString());
            }
        }



        /// <summary>
        /// メッセージをログに追加する。
        /// </summary>
        /// <param name="messages">メッセージ</param>
        public static void AddRange(string[] messages)
        {
            lock (_syncObj) _logs.AddRange(messages);
        }

        /// <summary>
        /// ログにセパレータを追加する。
        /// </summary>
        public static void AddSeparator()
        {
            lock (_syncObj) _logs.Add(_sepalator);
        }

        /// <summary>
        /// ログをクリアする。
        /// </summary>
        public static void Clear()
        {
            lock (_syncObj) _logs.Clear();
        }

        /// <summary>
        /// ログをファイルに保存する。
        /// </summary>
        /// <param name="logFilePath">ログファイルPath</param>
        public static void Save(string logFilePath)
        {
            using var writer = new StreamWriter(logFilePath);
            foreach (var log in _logs)
            {
                writer.WriteLine(log);
            }
        }
    }
}
