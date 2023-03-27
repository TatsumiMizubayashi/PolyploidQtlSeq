namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// 親2 SNP-indexルール
    /// </summary>
    internal class Parent2SnpIndexRule : IQtlSeqTargetVariantRule
    {
        private readonly Parent2SnpIndexRange _range;

        /// <summary>
        /// 親2 SNP-incexルールを作成する。
        /// </summary>
        /// <param name="range">範囲</param>
        public Parent2SnpIndexRule(Parent2SnpIndexRange range)
        {
            _range = range;
        }


        public bool Ok(SnpIndexVariant variant)
        {
            var snpIndex = variant.Parent2.SnpIndex;

            return _range.Lower <= snpIndex && snpIndex <= _range.Upper;
        }
    }
}
