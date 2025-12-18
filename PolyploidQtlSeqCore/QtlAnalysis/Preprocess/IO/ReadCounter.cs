namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// リード数カウンターの抽象クラス
    /// </summary>
    internal abstract class ReadCounter
    {
        /// <summary>
        /// リードカウンターを作成する。
        /// </summary>
        /// <param name="type">アレルタイプ</param>
        /// <returns>リードカウンター</returns>
        public static ReadCounter Create(VariantAlleleType type)
        {
            return type switch
            {
                VariantAlleleType.RefAlt => new RefAltReadCounter(),
                VariantAlleleType.Multi => new MultiAlleleReadCounter(),
                _ => throw new InvalidOperationException()
            };
        }

        /// <summary>
        /// リード数をカウントする。
        /// </summary>
        /// <param name="ad">AD</param>
        /// <returns>(refCount, altCount)</returns>
        public abstract (int refCount, int altCount) Count(AD ad);
    }
}
