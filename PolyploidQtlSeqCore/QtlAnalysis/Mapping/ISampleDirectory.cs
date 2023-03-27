using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QtlAnalysis.Mapping
{
    /// <summary>
    /// サンプルディレクトリインターフェイス
    /// </summary>
    public interface ISampleDirectory
    {
        /// <summary>
        /// ディレクトリのPathを取得する。
        /// </summary>
        public string Path { get; }

        /// <summary>
        /// サンプル名（ディレクトリ名）を取得する。
        /// </summary>
        public string SampleName { get; }

        /// <summary>
        /// Fastqファイルペア情報に変換する。
        /// </summary>
        /// <returns>FastqFilePair配列</returns>
        internal FastqFilePair[] ToFastqFilePairs()
        {
            return FastqFilePairEnumerator.Enumerate(Path);
        }

        /// <summary>
        /// ディレクトリ内にBAMファイルがあるかどうかを調査する。
        /// </summary>
        /// <returns>BAMファイルがある場合はtrue</returns>
        internal bool HasBamFile()
        {
            var bamFiles = BamFileEnumerator.Enumerate(Path);

            return bamFiles.Any();
        }

        /// <summary>
        /// BAMファイル情報に変換する。
        /// </summary>
        /// <returns>BAMファイル情報</returns>
        internal BamFile ToBamFile()
        {
            var bamFiles = BamFileEnumerator.Enumerate(Path);
            if (bamFiles.Length == 0) throw new InvalidOperationException("BAMファイルがありません。");
            if (bamFiles.Length > 1) throw new InvalidOperationException("BAMファイルが複数あります。");

            return new BamFile(SampleName, bamFiles[0]);
        }
    }
}
