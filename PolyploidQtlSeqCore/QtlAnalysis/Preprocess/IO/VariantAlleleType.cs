namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// 変異アレルの種類
    /// </summary>
    internal enum VariantAlleleType
    {
        /// <summary>
        /// RefアレルとAltアレルの組み合わせ
        /// </summary>
        RefAlt = 0,

        /// <summary>
        /// Alt1アレルとAlt2アレルの組み合わせ
        /// </summary>
        Alt1Alt2,

        /// <summary>
        /// ３種類以上の組み合わせ
        /// </summary>
        Multi
    }
}
