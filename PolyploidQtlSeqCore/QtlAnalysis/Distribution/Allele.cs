namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// アレル
    /// </summary>
    internal class Allele
    {
        /// <summary>
        /// アレルを作成する。
        /// </summary>
        /// <param name="isAlt">Alt型かどうか</param>
        public Allele(bool isAlt)
        {
            IsAlt = isAlt;
        }

        /// <summary>
        /// Alt型かどうかを取得する。
        /// </summary>
        public bool IsAlt { get; }
    }
}
