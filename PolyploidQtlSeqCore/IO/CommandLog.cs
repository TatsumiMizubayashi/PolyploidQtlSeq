namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// 実行したコマンドのログ
    /// </summary>
    internal static class CommandLog
    {
        private static readonly List<string> _commands = [];
        private static readonly object _syncObj = new();

        /// <summary>
        /// ログに実行したコマンドを追加する。
        /// </summary>
        /// <param name="command">コマンド</param>
        public static void Add(string command)
        {
            lock (_syncObj) _commands.Add(command);
        }

        /// <summary>
        /// コマンドログをクリアする。
        /// </summary>
        public static void Clear()
        {
            lock (_syncObj) _commands.Clear();
        }

        /// <summary>
        /// コマンドログをファイルに保存する。
        /// </summary>
        /// <param name="logFilePath">ログファイルPath</param>
        public static void Save(string logFilePath)
        {
            using var writer = new StreamWriter(logFilePath);
            foreach (var command in _commands)
            {
                writer.WriteLine(command);
                writer.WriteLine();
            }
        }
    }
}
