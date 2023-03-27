namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// Impact
    /// </summary>
    [Flags]
    internal enum Impact
    {
        /// <summary>
        /// Snp Eff情報なし
        /// </summary>
        None = 0,


        /// <summary>
        /// HIGH
        /// </summary>
        HIGH = 1 << 0,

        /// <summary>
        /// MODERATE
        /// </summary>
        MODERATE = 1 << 1,

        /// <summary>
        /// LOW
        /// </summary>
        LOW = 1 << 2,

        /// <summary>
        /// MODIFIER
        /// </summary>
        MODIFIER = 1 << 3,

    }
}
