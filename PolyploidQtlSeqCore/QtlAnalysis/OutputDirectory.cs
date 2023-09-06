namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL解析結果出力ディレクトリ
    /// </summary>
    internal class OutputDirectory
    {
        /// <summary>
        /// 出力ディレクトリを作成する。
        /// </summary>
        /// <param name="outputDirPath">出力ディレクトリPath</param>
        public OutputDirectory(string outputDirPath)
        {
            if (string.IsNullOrWhiteSpace(outputDirPath)) throw new ArgumentNullException(nameof(outputDirPath));

            Path = System.IO.Path.GetFullPath(outputDirPath);
        }

        /// <summary>
        /// 出力ディレクトリPATHを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// 出力ディレクトリ内のファイルPathを作成する。
        /// </summary>
        /// <param name="fileNames">ファイル名</param>
        /// <returns>出力ディレクトリ内ファイルPATH</returns>
        internal string CreateFilePath(params string[] fileNames)
        {
            var items = new[] { Path }.Concat(fileNames).ToArray();
            return System.IO.Path.Combine(items);
        }

        /// <summary>
        /// 出力ディレクトリを作成する。
        /// </summary>
        internal void Create()
        {
            if (Directory.Exists(Path)) return;

            Directory.CreateDirectory(Path);
        }
    }
}
