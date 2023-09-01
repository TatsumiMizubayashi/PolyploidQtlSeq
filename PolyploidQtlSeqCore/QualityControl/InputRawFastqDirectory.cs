using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 入力RawFastqディレクトリ
    /// </summary>
    internal class InputRawFastqDirectory
    {
        /// <summary>
        /// 入力RawFastqディレクトリを作成する。
        /// </summary>
        /// <param name="inputDirPath">入力RawFastqディレクトリPath</param>
        public InputRawFastqDirectory(string inputDirPath)
        {
            if (string.IsNullOrEmpty(inputDirPath)) throw new ArgumentException(null, nameof(inputDirPath));
            if (!Directory.Exists(inputDirPath)) throw new DirectoryNotFoundException($"{inputDirPath} not found.");

            Path = System.IO.Path.GetFullPath(inputDirPath);
        }

        /// <summary>
        /// 入力RawFastqディレクトリのPATHを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// ディレクトリ内のFastqファイルをFastqFilePair情報に変換する。
        /// </summary>
        /// <returns>Fastqファイルペア配列</returns>
        internal FastqFilePair[] ToFastqFilePairs()
        {
            return FastqFilePairEnumerator.Enumerate(Path);
        }
    }
}
