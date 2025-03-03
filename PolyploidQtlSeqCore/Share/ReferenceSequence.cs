namespace PolyploidQtlSeqCore.Share
{
    /// <summary>
    /// リファレンスシークエンス
    /// </summary>
    internal class ReferenceSequence
    {
        /// <summary>
        /// リファレンスシークエンスファイルを作成する。
        /// </summary>
        /// <param name="refSeqFilePath">リファレンスシークエンスファイルPath</param>
        public ReferenceSequence(string refSeqFilePath)
        {
            ArgumentException.ThrowIfNullOrEmpty(refSeqFilePath);
            if (!File.Exists(refSeqFilePath)) throw new FileNotFoundException($"{refSeqFilePath} not found.");

            Path = refSeqFilePath;
        }

         /// <summary>
        /// リファレンスシークエンスファイルPathを取得する。
        /// </summary>
        internal string Path { get; }


        /// <summary>
        /// bcftools用引数に変換する。
        /// </summary>
        /// <returns>bcf引数</returns>
        internal string ToBcftoolsArg()
        {
            return $"-f {Path}";
        }
    }
}
