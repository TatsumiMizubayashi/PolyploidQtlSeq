namespace Sequence.Position.Extract
{
    /// <summary>
    /// ゲノム位置を持つ項目のChunk.
    /// 同じ染色体名の連続した位置情報が複数まとめられている。
    /// </summary>
    /// <typeparam name="T">ゲノム位置情報を持つクラス</typeparam>
    internal sealed class GenomePositionChunk<T> : IHasGenomePositionItem
        where T : IHasGenomePositionItem
    {
        private readonly T[] _items;

        /// <summary>
        /// ゲノム位置を持つ項目のChunkを作成する。
        /// </summary>
        /// <param name="items">コレクション</param>
        public GenomePositionChunk(IEnumerable<T> items)
        {
            _items = [.. items];
            if (_items.Length == 0) throw new AggregateException($"{nameof(items)}が空です。");
            if (_items.Select(x => x.GenomePosition.ChrName).Distinct().Count() != 1)
                throw new ArgumentException($"異なる染色体名が存在します。");
            
            var areaStart = _items.Min(x => x.GenomePosition.Start);
            var areaEnd = _items.Max(x => x.GenomePosition.End);
            GenomePosition = new GenomePosition(_items[0].GenomePosition.ChrName, areaStart, areaEnd);
        }

        /// <summary>
        /// 所持しているゲノム位置を全て含む領域情報を取得する。
        /// </summary>
        public GenomePosition GenomePosition { get; }

        /// <summary>
        /// 指定した位置と位置が完全一致する項目を抽出する。
        /// </summary>
        /// <param name="targetPosition">位置</param>
        /// <returns>位置が完全一致する項目</returns>
        public T[] ExtractMatch(GenomePosition targetPosition)
        {
            return [.. _items.Where(x => x.GenomePosition.IsMatch(targetPosition))];
        }

        /// <summary>
        /// 指定した位置と重複する項目を抽出する。
        /// </summary>
        /// <param name="targetPosition">位置</param>
        /// <returns>位置が重複する項目</returns>
        public T[] ExtractOverlap(GenomePosition targetPosition)
        {
            return [.. _items.Where(x => x.GenomePosition.IsOverlap(targetPosition))];
        }
    }
}
